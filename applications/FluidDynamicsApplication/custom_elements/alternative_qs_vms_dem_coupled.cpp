//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

// Project includes
#include "includes/cfd_variables.h"
#include "includes/dem_variables.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

// Aplication includes
#include "alternative_qs_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/sparse_matrix_multiplication_utility.h"

namespace Kratos
{
DenseVector<std::vector<double>> mExactPorosity;
DenseVector<Matrix> mExactBodyForce;
DenseVector<Matrix> mExactPorosityGradient;
//////////////////////////Life cycle

template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId):
    QSVMS<TElementData>(NewId)
{}

template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    QSVMS<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    QSVMS<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    QSVMS<TElementData>(NewId,pGeometry,pProperties)
{}

///////////Destructor

template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::~AlternativeQSVMSDEMCoupled()
{}

template< class TElementData >
Element::Pointer AlternativeQSVMSDEMCoupled<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AlternativeQSVMSDEMCoupled>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< class TElementData >
Element::Pointer AlternativeQSVMSDEMCoupled<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AlternativeQSVMSDEMCoupled>(NewId, pGeom, pProperties);
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class does things with constitutive law here.
    QSVMS<TElementData>::Initialize(rCurrentProcessInfo);

    if(Dim == 2){
        if (NumNodes == 9 || NumNodes == 6)
            mInterpolationOrder = 2;
    }
    else if(Dim == 3){
        if (NumNodes == 10 || NumNodes == 27)
            mInterpolationOrder = 2;
    }

    const unsigned int number_of_gauss_points = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mPreviousVelocity.size() != number_of_gauss_points)
    {
        mPreviousVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mPreviousVelocity[g] = ZeroVector(Dim);
    }

    if (mPredictedSubscaleVelocity.size() != number_of_gauss_points)
    {
        mPredictedSubscaleVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mPredictedSubscaleVelocity[g] = ZeroVector(Dim);
    }
    // The prediction is updated before each non-linear iteration:
    // It is not stored in a restart and can be safely initialized.

    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mViscousResistanceTensor.size() != number_of_gauss_points)
    {
        mViscousResistanceTensor.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mViscousResistanceTensor[g] = ZeroMatrix(Dim,Dim);
    }
}


template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::GetShapeSecondDerivatives(
    DenseVector<DenseVector<Matrix>> &rDDN_DDX) const
{
    const GeometryData::IntegrationMethod integration_method = this->GetIntegrationMethod();
    const GeometryType& r_geometry = this->GetGeometry();

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = r_geometry.IntegrationPoints(integration_method);

    if (rDDN_DDX.size() != IntegrationPoints.size()){
        DenseVector<DenseVector<Matrix>> temp(IntegrationPoints.size());
        rDDN_DDX.swap(temp);
    }

    Matrix J(r_geometry.WorkingSpaceDimension(),r_geometry.LocalSpaceDimension());
    Matrix Jinv(r_geometry.LocalSpaceDimension(), r_geometry.WorkingSpaceDimension());
    DenseVector<Matrix> GradJ(r_geometry.LocalSpaceDimension());

    double DetJ;
    double DetA;
    const ShapeFunctionDerivativesArrayType DN_De = r_geometry.ShapeFunctionsLocalGradients(integration_method);
    for ( unsigned int g = 0; g < IntegrationPoints.size(); g++ ){

        DenseVector<Matrix> aux;

        if (aux.size() != r_geometry.PointsNumber()){
            DenseVector<Matrix> temp(r_geometry.PointsNumber());
            aux.swap( temp );
            rDDN_DDX[g].swap(temp);
        }

        Matrix DN_DX;
        if (DN_DX.size1() != r_geometry.PointsNumber() || DN_DX.size2() != r_geometry.LocalSpaceDimension())
            DN_DX.resize( r_geometry.PointsNumber(), r_geometry.LocalSpaceDimension(), false );


        const GeometryType::CoordinatesArrayType& local_point_coordinates = IntegrationPoints[g];

        ShapeFunctionsSecondDerivativesType DDN_DDe;
        r_geometry.ShapeFunctionsSecondDerivatives(DDN_DDe, local_point_coordinates);


        Matrix A, Ainv;

        r_geometry.Jacobian(J,g,integration_method);
        MathUtils<double>::InvertMatrix( J, Jinv, DetJ );

        DN_DX = prod(DN_De[g],Jinv);

        if(Dim == 2){
            A = ZeroMatrix(3,3);
            Ainv = ZeroMatrix(3,3);

            A(0,0) = J(0,0) * J(0,0);
            A(0,1) = J(1,0) * J(1,0);
            A(0,2) = 2.0 * J(0,0) * J(1,0);

            A(1,0) = J(0,1) * J(0,1);
            A(1,1) = J(1,1) * J(1,1);
            A(1,2) = 2.0 * J(0,1) * J(1,1);

            A(2,0) = J(0,0) * J(0,1);
            A(2,1) = J(1,0) * J(1,1);
            A(2,2) = J(0,0) * J(1,1) + J(0,1) * J(1,0);
        }
        else if(Dim == 3){
            A = ZeroMatrix(6,6);
            Ainv = ZeroMatrix(6,6);

            A(0,0) = J(0,0) * J(0,0);
            A(0,1) = J(1,0) * J(1,0);
            A(0,2) = J(2,0) * J(2,0);
            A(0,3) = 2.0 * J(0,0) * J(1,0);
            A(0,4) = 2.0 * J(1,0) * J(2,0);
            A(0,5) = 2.0 * J(0,0) * J(2,0);

            A(1,0) = J(0,1) * J(0,1);
            A(1,1) = J(1,1) * J(1,1);
            A(1,2) = J(2,1) * J(2,1);
            A(1,3) = 2.0 * J(0,1) * J(1,1);
            A(1,4) = 2.0 * J(1,1) * J(2,1);
            A(1,5) = 2.0 * J(0,1) * J(2,1);

            A(2,0) = J(0,2) * J(0,2);
            A(2,1) = J(1,2) * J(1,2);
            A(2,2) = J(2,2) * J(2,2);
            A(2,3) = 2.0 * J(0,2) * J(1,2);
            A(2,4) = 2.0 * J(1,2) * J(2,2);
            A(2,5) = 2.0 * J(0,2) * J(2,2);

            A(3,0) = J(0,0) * J(0,1);
            A(3,1) = J(1,0) * J(1,1);
            A(3,2) = J(2,0) * J(2,1);
            A(3,3) = J(0,0) * J(1,1) + J(0,1) * J(1,0);
            A(3,4) = J(1,0) * J(2,1) + J(1,1) * J(2,0);
            A(3,5) = J(0,0) * J(2,1) + J(0,1) * J(2,0);

            A(4,0) = J(0,1) * J(0,2);
            A(4,1) = J(1,1) * J(1,2);
            A(4,2) = J(2,1) * J(2,2);
            A(4,3) = J(0,1) * J(1,2) + J(0,2) * J(1,1);
            A(4,4) = J(1,1) * J(2,2) + J(1,2) * J(2,1);
            A(4,5) = J(0,1) * J(2,2) + J(0,2) * J(2,1);

            A(5,0) = J(0,0) * J(0,2);
            A(5,1) = J(1,0) * J(1,2);
            A(5,2) = J(2,0) * J(2,2);
            A(5,3) = J(0,0) * J(1,2) + J(0,2) * J(1,0);
            A(5,4) = J(1,0) * J(2,2) + J(1,2) * J(2,0);
            A(5,5) = J(0,0) * J(2,2) + J(0,2) * J(2,0);

        }

        MathUtils<double>::InvertMatrix( A, Ainv, DetA );
        DenseVector<Matrix> H(r_geometry.WorkingSpaceDimension());
        for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); d++)
            H[d] = ZeroMatrix(r_geometry.LocalSpaceDimension(),r_geometry.LocalSpaceDimension());

        for (IndexType p = 0; p < r_geometry.PointsNumber(); ++p) {
            const array_1d<double, 3>& r_coordinates = r_geometry[p].Coordinates();
            H[0](0,0) += r_coordinates[0] * DDN_DDe[p](0,0);
            H[0](0,1) += r_coordinates[0] * DDN_DDe[p](0,1);
            H[0](1,0) += r_coordinates[0] * DDN_DDe[p](1,0);
            H[0](1,1) += r_coordinates[0] * DDN_DDe[p](1,1);

            H[1](0,0) += r_coordinates[1] * DDN_DDe[p](0,0);
            H[1](0,1) += r_coordinates[1] * DDN_DDe[p](0,1);
            H[1](1,0) += r_coordinates[1] * DDN_DDe[p](1,0);
            H[1](1,1) += r_coordinates[1] * DDN_DDe[p](1,1);

            if(Dim == 3){
                H[0](0,2) += r_coordinates[0] * DDN_DDe[p](0,2);
                H[0](1,2) += r_coordinates[0] * DDN_DDe[p](1,2);
                H[0](2,0) += r_coordinates[0] * DDN_DDe[p](2,0);
                H[0](2,1) += r_coordinates[0] * DDN_DDe[p](2,1);
                H[0](2,2) += r_coordinates[0] * DDN_DDe[p](2,2);

                H[1](0,2) += r_coordinates[1] * DDN_DDe[p](0,2);
                H[1](1,2) += r_coordinates[1] * DDN_DDe[p](1,2);
                H[1](2,0) += r_coordinates[1] * DDN_DDe[p](2,0);
                H[1](2,1) += r_coordinates[1] * DDN_DDe[p](2,1);
                H[1](2,2) += r_coordinates[1] * DDN_DDe[p](2,2);

                H[2](0,0) += r_coordinates[2] * DDN_DDe[p](0,0);
                H[2](0,1) += r_coordinates[2] * DDN_DDe[p](0,1);
                H[2](1,0) += r_coordinates[2] * DDN_DDe[p](1,0);
                H[2](1,1) += r_coordinates[2] * DDN_DDe[p](1,1);
                H[2](0,2) += r_coordinates[2] * DDN_DDe[p](0,2);
                H[2](1,2) += r_coordinates[2] * DDN_DDe[p](1,2);
                H[2](2,0) += r_coordinates[2] * DDN_DDe[p](2,0);
                H[2](2,1) += r_coordinates[2] * DDN_DDe[p](2,1);
                H[2](2,2) += r_coordinates[2] * DDN_DDe[p](2,2);
            }
        }

        for (IndexType p = 0; p < r_geometry.PointsNumber(); ++p) {
            Vector result, rhs;
            if (Dim == 2){
                rhs = ZeroVector(3);
                result = ZeroVector(3);
                rhs[0] = DDN_DDe[p](0,0) - DN_DX(p,0) * H[0](0,0) - DN_DX(p,1) * H[1](0,0);
                rhs[1] = DDN_DDe[p](1,1) - DN_DX(p,0) * H[0](1,1) - DN_DX(p,1) * H[1](1,1);
                rhs[2] = DDN_DDe[p](0,1) - DN_DX(p,0) * H[0](0,1) - DN_DX(p,1) * H[1](0,1);
            }
            else if (Dim == 3){
                rhs = ZeroVector(6);
                result = ZeroVector(6);
                rhs[0] = DDN_DDe[p](0,0) - DN_DX(p,0) * H[0](0,0) - DN_DX(p,1) * H[1](0,0) - DN_DX(p,2) * H[2](0,0);
                rhs[1] = DDN_DDe[p](1,1) - DN_DX(p,0) * H[0](1,1) - DN_DX(p,1) * H[1](1,1) - DN_DX(p,2) * H[2](1,1);
                rhs[2] = DDN_DDe[p](2,2) - DN_DX(p,0) * H[0](2,2) - DN_DX(p,1) * H[1](2,2) - DN_DX(p,2) * H[2](2,2);
                rhs[3] = DDN_DDe[p](0,1) - DN_DX(p,0) * H[0](0,1) - DN_DX(p,1) * H[1](0,1) - DN_DX(p,2) * H[2](0,1);
                rhs[4] = DDN_DDe[p](1,2) - DN_DX(p,0) * H[0](1,2) - DN_DX(p,1) * H[1](1,2) - DN_DX(p,2) * H[2](1,2);
                rhs[5] = DDN_DDe[p](0,2) - DN_DX(p,0) * H[0](0,2) - DN_DX(p,1) * H[1](0,2) - DN_DX(p,2) * H[2](0,2);
            }

            aux[p].resize(r_geometry.WorkingSpaceDimension(), r_geometry.WorkingSpaceDimension(), false );

            result = prod(Ainv, rhs);

            aux[p](0,0) = result[0];
            aux[p](1,1) = result[1];
            if (Dim == 2){
                aux[p](0,1) = result[2];
                aux[p](1,0) = result[2];
            }
            else if (Dim == 3){
                aux[p](2,2) = result[2];
                aux[p](0,1) = result[3];
                aux[p](1,0) = result[3];
                aux[p](0,2) = result[5];
                aux[p](2,0) = result[5];
                aux[p](2,1) = result[4];
                aux[p](1,2) = result[4];
            }
        }

        rDDN_DDX[g] = aux;
    }
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    const SizeType number_of_integration_points = integration_points.size();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(
        gauss_weights, shape_functions, shape_derivatives);

    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);
    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    for (unsigned int g = 0; g < number_of_integration_points; g++ ) {
        data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),shape_derivatives[g]);
        if (rVariable == PRESSURE) {
            const auto& r_pressure = data.Pressure;
            double value = this->GetAtCoordinate(r_pressure,data.N);
            rOutput[g] = value;
        }
    }
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    const SizeType number_of_integration_points = integration_points.size();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(
        gauss_weights, shape_functions, shape_derivatives);

    if (rValues.size() != number_of_integration_points)
        rValues.resize(number_of_integration_points);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++ ) {
        data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),shape_derivatives[g]);
        Matrix value = ZeroMatrix(Dim, Dim);
        if (rVariable == VELOCITY_GRADIENT) {
            const auto& r_velocity = data.Velocity;
            for (unsigned int i = 0; i < NumNodes; i++) {
                for (unsigned int d = 0; d < Dim; d++) {
                    for (unsigned int e = 0; e < Dim; e++)
                        value(d,e) += data.DN_DX(i,d) * r_velocity(i,e);
                }
            }
        }
        rValues[g] = value;
    }
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3>>& rVariable,
    std::vector<array_1d<double,3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    const SizeType number_of_integration_points = integration_points.size();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(
        gauss_weights, shape_functions, shape_derivatives);

    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++ ) {
        data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),shape_derivatives[g]);
        array_1d<double,3> value(3,0.0);
        if (rVariable == VELOCITY) {
            const auto& r_velocity = data.Velocity;
            value = this->GetAtCoordinate(r_velocity,data.N);
        }
        if (rVariable == BODY_FORCE) {
            value = this->GetAtCoordinate(data.BodyForce,data.N);
        }
        if (rVariable == PRESSURE_GRADIENT){
            const auto& r_pressure = data.Pressure;
            for (unsigned int i = 0; i < NumNodes; i++) {
                for (unsigned int d = 0; d < Dim; d++) {
                    value[d] += r_pressure[i] * data.DN_DX(i,d);
                }
            }
        }
        rOutput[g] = value;
    }
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    }
    else if (rVariable == SUBSCALE_VELOCITY){
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        this->GetShapeSecondDerivatives(shape_function_second_derivatives);
        const unsigned int NumGauss = GaussWeights.size();
        array_1d<double,NumNodes*Dim> momentum_rhs = ZeroVector(NumNodes*Dim);
        VectorType MassRHS = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        for (unsigned int g = 0; g < NumGauss; g++)
        {

            this->UpdateIntegrationPointData(data, g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g], shape_function_second_derivatives[g]);

            array_1d<double, 3> MomentumRes = ZeroVector(3);
            double MassRes = 0.0;

            //array_1d<double,3> convective_velocity = this->GetAtCoordinate(data.Velocity,data.N) - this->GetAtCoordinate(data.MeshVelocity,data.N);
            array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(data);
            this->MomentumProjTerm(data, convective_velocity, MomentumRes);
            this->MassProjTerm(data,MassRes);

        //     for (unsigned int i = 0; i < NumNodes; i++)
        //     {
        //         array_1d<double,3>& momentum_rhs = this->GetGeometry()[i].GetValue(ADVPROJ);
        //         VectorType MassRHS = this->GetGeometry()[i].GetValue(DIVPROJ);

        //         double W = data.Weight*data.N[i];
        //         unsigned int row = i*Dim;
        //         for (unsigned int d = 0; d < Dim; d++)
        //             momentum_rhs[row+d] += W*MomentumRes[d];
        //         NodalArea[i] += W;
        //         MassRHS[i] += W*MassRes;
        //         }
        // }
        /* Projections of the elemental residual are computed with
            * Newton-Raphson iterations of type M(lumped) dx = ElemRes - M(consistent) * x
        */

        // Carefully write results to nodal variables, to avoid parallelism problems
            for (unsigned int i = 0; i < NumNodes; ++i)
            {
                double W = data.Weight*data.N[i];
                // Write nodal area
                unsigned int row = i*Dim;
                //auto& momentum_rhs = this->GetGeometry()[i].GetValue(ADVPROJ);
                //double& MassRHS = this->GetGeometry()[i].GetValue(DIVPROJ);
                for (unsigned int d = 0; d < Dim; d++)
                    momentum_rhs[row+d] += data.N[i]*MomentumRes[d];
                NodalArea[i] += W;
                MassRHS[i] += data.N[i]*MassRes;
                }
        }

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            unsigned int row_i = i*Dim;
            double W = data.Weight*data.N[i];
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            auto& rMomValue = this->GetGeometry()[i].GetValue(ADVPROJ);
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            // Substract M(consistent)*x(i-1) from RHS
            for(unsigned int j = 0; j < NumNodes; ++j) // RHS -= Weigth * Ones(TNumNodes,TNumNodes) * x(i-1)
            {
                unsigned int row_j = j*Dim;
                for(unsigned int d = 0; d < Dim; ++d)
                    momentum_rhs[row_j+d] -= W * this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ)[d];
                MassRHS[j] -= W * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
            }
            for(unsigned int d = 0; d < Dim; ++d){ // RHS -= Weigth * Identity(TNumNodes,TNumNodes) * x(i-1)
                momentum_rhs[row_i+d] -= W * this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ)[d];
                rMomValue[d] += momentum_rhs[row_i+d];
            }
            MassRHS[i] -= W * this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
            this->GetGeometry()[i].GetValue(DIVPROJ) += MassRHS[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        }
    }
}

template< class TElementData >
GeometryData::IntegrationMethod AlternativeQSVMSDEMCoupled<TElementData>::GetIntegrationMethod() const
{
    if(mInterpolationOrder == 1)
        return GeometryData::IntegrationMethod::GI_GAUSS_2;
    else
        return GeometryData::IntegrationMethod::GI_GAUSS_3;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int AlternativeQSVMSDEMCoupled<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    int out = QSVMS<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    const auto &r_geom = this->GetGeometry();
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const auto& rNode = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);
    }

    return out;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
std::string AlternativeQSVMSDEMCoupled<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "AlternativeQSVMSDEMCoupled #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "AlternativeQSVMSDEMCoupled" << Dim << "D";
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();
    this->GetShapeSecondDerivatives(shape_function_second_derivatives);
    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointData(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g],shape_function_second_derivatives[g]);

        this->CalculateResistanceTensor(data);
        //this->UpdateSubscaleVelocity(data);
    }
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    //Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();
    this->GetShapeSecondDerivatives(shape_function_second_derivatives);
    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointData(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g],shape_function_second_derivatives[g]);

        this->UpdateSubscaleVelocity(data);
    }
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX,
    const typename TElementData::ShapeFunctionsSecondDerivativesType& rDDN_DDX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
    rData.UpdateSecondDerivativesValues(rDDN_DDX);

    this->CalculateMaterialResponse(rData);
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3> &rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    const auto& rGeom = this->GetGeometry();

    Vector convection; // u * grad(N)
    this->ConvectionOperator(convection,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];
    const auto& body_force = this->GetAtCoordinate(rData.BodyForce, rData.N);

    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;
    array_1d<double,Dim> fluid_fraction_gradient = ZeroVector(Dim);
    for (unsigned int i = 0; i < NumNodes; i++)
        for (unsigned int d = 0; d < Dim; d++)
            fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];

    Vector grad_alpha_sym_grad_u, sigma_U, grad_div_u, div_sym_grad_u;
    BoundedMatrix<double,Dim,Dim> sym_gradient_u;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        sigma_U = ZeroVector(Dim);
        grad_div_u = ZeroVector(Dim);
        sym_gradient_u = ZeroMatrix(Dim, Dim);
        grad_alpha_sym_grad_u = ZeroVector(Dim);
        div_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            double div_u = 0.0;
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
                sym_gradient_u(d,e) += 1.0/2.0 * (rData.DN_DX(i,d) * r_velocities(i,e) + rData.DN_DX(i,e) * r_velocities(i,d));
                grad_alpha_sym_grad_u[d] += fluid_fraction_gradient[e] * sym_gradient_u(d,e);
                div_u += rData.DN_DX(i,e) * r_velocities(i,e);
                grad_div_u[d] += rData.DDN_DDX[i](d,e) *  r_velocities(i,e);
                if (d == e)
                    div_sym_grad_u[d] += rData.DDN_DDX[i](e,e) * r_velocities(i,d);
                else
                    div_sym_grad_u[d] += 1.0/2.0 * (rData.DDN_DDX[i](e,d) * r_velocities(i,e) + rData.DDN_DDX[i](e,e) * r_velocities(i,d));
            }
            rResidual[d] += density * (/*rData.N[i] * r_body_forces(i,d)*/ - fluid_fraction * rData.N[i] * r_acceleration[d] - fluid_fraction * convection[i] * r_velocities(i,d)) + 2.0 * grad_alpha_sym_grad_u[d] * viscosity - 2.0 / 3.0 * viscosity * fluid_fraction_gradient[d] * div_u + 2.0 * fluid_fraction * viscosity * div_sym_grad_u[d] - 2.0/3.0 * fluid_fraction * viscosity * grad_div_u[d] - fluid_fraction * rData.DN_DX(i,d) * r_pressures[i] - sigma_U[d];
            }
        }
    for (unsigned int d = 0; d < Dim; d++)
        rResidual[d] += density * body_force[d];
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::MomentumProjTerm(
    const TElementData& rData,
    const array_1d<double,3>& rConvectionVelocity,
    array_1d<double,3> &rMomentumRHS) const
{
    //const auto& rGeom = this->GetGeometry();

    Vector AGradN;
    this->ConvectionOperator(AGradN,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const auto body_force = this->GetAtCoordinate(rData.BodyForce, rData.N);

    const auto r_velocities = rData.Velocity;
    const auto r_pressures = rData.Pressure;
    array_1d<double,3> fluid_fraction_gradient = ZeroVector(Dim);
    for (unsigned int i = 0; i < NumNodes; i++)
        for (unsigned int d = 0; d < Dim; d++)
            fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];

    Vector grad_alpha_sym_grad_u, grad_div_u, div_sym_grad_u;
    BoundedMatrix<double,Dim,Dim> sym_gradient_u;
    for (unsigned int i = 0; i < NumNodes; i++) {
        grad_div_u = ZeroVector(Dim);
        sym_gradient_u = ZeroMatrix(Dim, Dim);
        grad_alpha_sym_grad_u = ZeroVector(Dim);
        div_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            double div_u = 0.0;
            for (unsigned int e = 0; e < Dim; e++){
                sym_gradient_u(d,e) += 1.0/2.0 * (rData.DN_DX(i,d) * r_velocities(i,e) + rData.DN_DX(i,e) * r_velocities(i,d));
                grad_alpha_sym_grad_u[d] += fluid_fraction_gradient[e] * sym_gradient_u(d,e);
                div_u += rData.DN_DX(i,e) * r_velocities(i,e);
                grad_div_u[d] += rData.DDN_DDX[i](d,e) * r_velocities(i,e);
                if (d == e)
                    div_sym_grad_u[d] += rData.DDN_DDX[i](e,e) * r_velocities(i,d);
                else
                    div_sym_grad_u[d] += 1.0/2.0 * (rData.DDN_DDX[i](e,d) * r_velocities(i,e) + rData.DDN_DDX[i](e,e) * r_velocities(i,d));
            }
            rMomentumRHS[d] += density * (- fluid_fraction * AGradN[i] * r_velocities(i,d)) + 2.0 * grad_alpha_sym_grad_u[d] * viscosity - 2.0/3.0 * viscosity * fluid_fraction_gradient[d] * div_u + 2.0 * fluid_fraction * viscosity * div_sym_grad_u[d] - 2.0/3.0 * fluid_fraction * viscosity * grad_div_u[d] - fluid_fraction * rData.DN_DX(i,d) * r_pressures[i];
        }
    }
    for (unsigned int d = 0; d < Dim; d++)
        rMomentumRHS[d] += density * body_force[d];
}

template<class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType& rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);

    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying convective operator by density to have correct units
    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    array_1d<double,3> fluid_fraction_gradient = ZeroVector(Dim);
    for (unsigned int i = 0; i < NumNodes; i++)
        for (unsigned int d = 0; d < Dim; d++)
            fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    double W = rData.Weight * density; // This density is for the dynamic term in the residual (rho*Du)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++)
            {
                // grad(q) * TauOne * du
                double UGAlpha = tau_one(d,d) * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                // u*grad(v) * TauOne * du
                double AU = tau_one(d,d) * fluid_fraction * fluid_fraction * AGradN[i] * rData.N[j];

                for (unsigned int e = 0; e < Dim; ++e){
                    double LI = tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](d,e) * rData.N[j];
                    double CI = 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](d,e) * rData.N[j];
                    double RSigmaU = tau_one(d,d) * fluid_fraction * sigma(d,e) * rData.N[i] * rData.N[j];
                    double DBetaU = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * rData.DN_DX(i,d) * rData.N[j] * fluid_fraction_gradient[e];
                    double GBetaU = tau_one(d,d) * fluid_fraction * viscosity * fluid_fraction_gradient[d] * rData.DN_DX(i,e) * rData.N[j];
                    for (unsigned int f = 0; f < Dim; ++f){
                        if (d == e){
                            LI += tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](f,f) * rData.N[j];
                            GBetaU += tau_one(d,d) * fluid_fraction * viscosity * fluid_fraction_gradient[f] * rData.DN_DX(i,f) * rData.N[j];
                        }
                    }

                    rMassMatrix(row+d, col+e) += W * (GBetaU + LI - CI - RSigmaU - DBetaU);
                }
                rMassMatrix(row+d,col+d) += W * (AU);
                rMassMatrix(row+Dim,col+d) += W * UGAlpha;
            }
        }
    }
}

template<class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::AddReactionStabilization(
    TElementData& rData,
    auto& rLHS,
    VectorType& rLocalRHS)
{

    const double density = this->GetAtCoordinate(rData.Density, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce, rData.N);

    Vector AGradN;
    this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX); // Get a * grad(Ni)

    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    array_1d<double,3> fluid_fraction_gradient = ZeroVector(Dim);
    for (unsigned int i = 0; i < NumNodes; i++)
        for (unsigned int d = 0; d < Dim; d++)
            fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];

    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);

    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        // Loop over columns
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++) // iterate over dimensions for velocity Dofs in this node combination
            {
                double RSigmaG = 0.0;
                double GR = 0.0;
                for (unsigned int e = 0; e < Dim; e++){
                    double ASigma = tau_one(d,d) * fluid_fraction * AGradN[i] * rData.N[j] * sigma(d,e);
                    double RSigmaA = tau_one(d,d) * fluid_fraction * rData.N[i] * AGradN[j] * sigma(d,e);
                    double LSigma_1 = 0.0;
                    double LSigma_2 = 0.0;
                    double CSigma = 0.0;
                    double GBetaSigma_1 = 0.0;
                    double GBetaSigma_2 = 0.0;
                    double DBetaSigma = 0.0;
                    double RSigmaL_1 = 0.0;
                    double RSigmaL_2 = 0.0;
                    double RSigmaC = 0.0;
                    double RGBeta_1 = 0.0;
                    double RGBeta_2 = 0.0;
                    double RDBeta = 0.0;
                    double RRSigma = 0.0;
                    GR += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,e) * sigma(e,d) * rData.N[j];
                    for (unsigned int f = 0; f < Dim; f++){
                        LSigma_1 += tau_one(d,d) * fluid_fraction * viscosity * rData.DDN_DDX[i](f,f) * sigma(d,e) * rData.N[j];
                        LSigma_2 += tau_one(d,d) * fluid_fraction * viscosity * rData.DDN_DDX[i](d,f) * sigma(f,e) * rData.N[j];
                        CSigma += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * rData.DDN_DDX[i](f,d) * sigma(f,e) * rData.N[j];
                        GBetaSigma_1 += tau_one(d,d) * viscosity * fluid_fraction_gradient[d] * rData.DN_DX(i,f) * sigma(f,e) * rData.N[j];
                        GBetaSigma_2 += tau_one(d,d) * viscosity * fluid_fraction_gradient[f] * rData.DN_DX(i,f) * sigma(d,e) * rData.N[j];
                        RSigmaL_1 += tau_one(d,d) * fluid_fraction * viscosity * rData.N[i] * sigma(d,e) * rData.DDN_DDX[j](f,f);
                        RSigmaL_2 += tau_one(d,d) * fluid_fraction * viscosity * rData.N[i] * sigma(d,f) * rData.DDN_DDX[j](e,f);
                        RSigmaC += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * sigma(d,f) * rData.N[i] * rData.DDN_DDX[j](f,e);
                        RGBeta_1 += tau_one(d,d) * viscosity * rData.N[i] * fluid_fraction_gradient[e] * sigma(d,f) * rData.DN_DX(j,f);
                        RGBeta_2 += tau_one(d,d) * viscosity * rData.N[i] * fluid_fraction_gradient[f] * rData.DN_DX(j,f) * sigma(d,e);
                        RDBeta += 2.0/3.0 * tau_one(d,d) * viscosity * rData.N[i] * sigma(d,f) * fluid_fraction_gradient[f] * rData.DN_DX(j,e);
                        RRSigma += tau_one(d,d) * sigma(d,f) * rData.N[i] * sigma(f,e) * rData.N[j];
                        DBetaSigma += 2.0 / 3.0 * tau_one(d,d) * viscosity * rData.DN_DX(i,d) * rData.N[j] * fluid_fraction_gradient[f] * sigma(f,e);
                    }
                    double LSigma = LSigma_1 + LSigma_2;
                    double GBetaSigma = GBetaSigma_1 + GBetaSigma_2;
                    double RSigmaL = RSigmaL_1 + RSigmaL_2;
                    double RGBeta = RGBeta_1 + RGBeta_2;
                    RSigmaG += tau_one(d,d) * fluid_fraction * sigma(d,e) * rData.N[i] * rData.DN_DX(j,e);
                    rLHS(row+d,col+e) += rData.Weight * (GBetaSigma + RGBeta - DBetaSigma - RDBeta + ASigma - RRSigma - RSigmaA + LSigma - CSigma + RSigmaL - RSigmaC);
                }
                rLHS(row+Dim,col+d) += rData.Weight * (GR);
                rLHS(row+d,col+Dim) += rData.Weight * (-RSigmaG);
            }
        }

        // RHS terms
        for (unsigned int d = 0; d < Dim; ++d)
        {
            double RSigmaF = 0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                RSigmaF += tau_one(d,d) * rData.N[i] * sigma(d,e) * body_force[e]; /*- momentum_projection[e]*/ //momentum_projection 0 because is ASGS
            }
            rLocalRHS[row+d] += rData.Weight * (- RSigmaF);
        }
    }
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::AddViscousTerm(
    const TElementData& rData,
    BoundedMatrix<double,LocalSize,LocalSize>& rLHS,
    VectorType& rRHS) {

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    BoundedMatrix<double,StrainSize,LocalSize> strain_matrix = ZeroMatrix(StrainSize,LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX,strain_matrix);

    const auto& constitutive_matrix = rData.C;
    BoundedMatrix<double,StrainSize,LocalSize> shear_stress_matrix = prod(constitutive_matrix,strain_matrix);

    // Multiply times integration point weight (I do this here to avoid a temporal in LHS += weight * Bt * C * B)
    strain_matrix *= rData.Weight;

    noalias(rLHS) += prod(trans(strain_matrix), fluid_fraction * shear_stress_matrix);
    noalias(rRHS) -= prod(trans(strain_matrix), fluid_fraction * rData.ShearStress);
}

// Add a the contribution from a single integration point to the velocity contribution
template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();
    //const auto& rGeom = this->GetGeometry();
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce,rData.N); // Force per unit of volume


    // const array_1d<double,3> convective_velocity =
    //     this->GetAtCoordinate(rData.Velocity, rData.N) -
    //     this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    const array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);
    array_1d<double,3> velocity = this->GetAtCoordinate(rData.Velocity,rData.N);

    array_1d<double,Dim>& r_prev_velocity = mPreviousVelocity[rData.IntegrationPointIndex];

    for (unsigned int n = 0; n < Dim; n++)
        r_prev_velocity[n] = velocity[n];

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // These two should be zero unless we are using OSS
    array_1d<double,3> MomentumProj = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    double MassProj = this->GetAtCoordinate(rData.MassProjection,rData.N);

    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    array_1d<double,Dim> fluid_fraction_gradient = ZeroVector(Dim);

    for (unsigned int i = 0; i < NumNodes; i++)
        for (unsigned int d = 0; d < Dim; d++)
            fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];

    AGradN *= density; // Convective term is always multiplied by density

    // Multiplying convective operator by density to have correct units
    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {

        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j

            // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            double V = fluid_fraction * rData.N[i] * AGradN[j];

            // q-p stabilization block (reset result)
            double G = 0.0;
            for (unsigned int d = 0; d < Dim; d++) {
                // Stabilization: u*grad(v) * TauOne * u*grad(u)
                // The last term comes from vh*d(u_ss)
                double AA = tau_one(d,d) * AGradN[i] * std::pow(fluid_fraction, 2) * AGradN[j];

                LHS(row+d,col+d) += rData.Weight * (V + AA);
                // Galerkin pressure term: Div(v) * p
                double P = fluid_fraction * rData.DN_DX(i,d) * rData.N[j];

                double QD = fluid_fraction * rData.DN_DX(j,d) * rData.N[i];

                double GP = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];

                double GAlphaD = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                G += tau_one(d,d) * std::pow(fluid_fraction,2) * rData.DN_DX(i,d) * rData.DN_DX(j,d);
                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                double GG_1 = 0.0;
                double GG_2 = 0.0;
                double DG = 0.0;
                double GGBeta_1 = 0.0;
                double GGBeta_2 = 0.0;
                double GC = 0.0;
                double CG = 0.0;
                double GL = 0.0;
                double LG = 0.0;
                // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                double AG = tau_one(d,d) * std::pow(fluid_fraction,2) * AGradN[i] * rData.DN_DX(j,d);
                double GA = tau_one(d,d) * std::pow(fluid_fraction,2) * rData.DN_DX(i,d) * AGradN[j];
                double GDBeta = 0.0;
                for (unsigned int e = 0; e < Dim; e++){
                    double DnuD = 2.0/3.0 * fluid_fraction * viscosity * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double GS = fluid_fraction * viscosity * rData.DN_DX(i,e) * rData.DN_DX(j,d);
                    GDBeta += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * fluid_fraction_gradient[e] * rData.DN_DX(i,e) * rData.DN_DX(j,d);
                    GL += tau_one(d,d) * viscosity * std::pow(fluid_fraction,2) * rData.DN_DX(i,e) * rData.DDN_DDX[j](d,e);
                    GL += tau_one(d,d) * viscosity * std::pow(fluid_fraction,2) * rData.DN_DX(i,d) * rData.DDN_DDX[j](e,e);
                    GGBeta_1 += tau_one(d,d) * fluid_fraction * viscosity * rData.DN_DX(i,e) * rData.DN_DX(j,e) * fluid_fraction_gradient[d];
                    GGBeta_2 += tau_one(d,d) * fluid_fraction * viscosity * rData.DN_DX(i,d) * rData.DN_DX(j,e) * fluid_fraction_gradient[e];
                    GG_1 += tau_one(d,d) * fluid_fraction * viscosity * fluid_fraction_gradient[d] * rData.DN_DX(i,e) * rData.DN_DX(j,e);
                    GG_2 += tau_one(d,d) * fluid_fraction * viscosity * fluid_fraction_gradient[e] * rData.DN_DX(i,e) * rData.DN_DX(j,d);
                    LG += tau_one(d,d) * viscosity * std::pow(fluid_fraction,2) * rData.DN_DX(j,e) * rData.DDN_DDX[i](d,e);
                    LG += tau_one(d,d) * viscosity * std::pow(fluid_fraction,2) * rData.DN_DX(j,d) * rData.DDN_DDX[i](e,e);
                    GC += 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DN_DX(i,e) * rData.DDN_DDX[j](e,d);
                    CG += 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](d,e) * rData.DN_DX(j,e);
                    DG += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.DN_DX(j,e);
                    double AL = tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * AGradN[i] * rData.DDN_DDX[j](d,e);
                    double LA = tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](d,e) * AGradN[j];
                    double AC = 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * AGradN[i] * rData.DDN_DDX[j](d,e);
                    double CA = 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](d,e) * AGradN[j];
                    double DBetaA = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * AGradN[j] * fluid_fraction_gradient[e] * rData.DN_DX(i,d);
                    double RSigma = rData.N[i] * sigma(d,e) * rData.N[j];
                    double AGBeta = tau_one(d,d) * fluid_fraction * viscosity * AGradN[i] * rData.DN_DX(j,d) * fluid_fraction_gradient[e];
                    double ADBeta = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * AGradN[i] * fluid_fraction_gradient[d] * rData.DN_DX(j,e);
                    double DD = tau_two * std::pow(fluid_fraction,2) * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double DU = tau_two * fluid_fraction * fluid_fraction_gradient[e] * rData.DN_DX(i,d) * rData.N[j];
                    double GD = tau_two * fluid_fraction * fluid_fraction_gradient[d] * rData.DN_DX(j,e) * rData.N[i];
                    double GU = tau_two * rData.N[j] * rData.N[i] * fluid_fraction_gradient[d] * fluid_fraction_gradient[e];
                    double GBetaA = tau_one(d,d) * fluid_fraction * viscosity * AGradN[j] * fluid_fraction_gradient[d] * rData.DN_DX(i,e);
                    double LL_diag_1 = 0.0;
                    double LL_diag_2 = 0.0;
                    double LL_2 = 0.0;
                    double LL_3 = 0.0;
                    double LL_4 = 0.0;
                    double LC_1 = 0.0;
                    double LC_2 = 0.0;
                    double LGBeta_1 = 0.0;
                    double LGBeta_2 = 0.0;
                    double LGBeta_3 = 0.0;
                    double LGBeta_4 = 0.0;
                    double LGBeta_5 = 0.0;
                    double LDBeta_1 = 0.0;
                    double LDBeta_2 = 0.0;
                    double CL_1 = 0.0;
                    double CL_2 = 0.0;
                    double CC = 0.0;
                    double CGBeta_1 = 0.0;
                    double CGBeta_2 = 0.0;
                    double CDBeta = 0.0;
                    double GBetaL_1 = 0.0;
                    double GBetaL_2 = 0.0;
                    double GBetaL_3 = 0.0;
                    double GBetaL_4 = 0.0;
                    double GBetaL_5 = 0.0;
                    double GBetaC_1 = 0.0;
                    double GBetaC_2 = 0.0;
                    double GBetaG_1 = 0.0;
                    double GBetaG_2 = 0.0;
                    double GBetaG_3 = 0.0;
                    double GBetaG_4 = 0.0;
                    double GBetaG_5 = 0.0;
                    double GBetaD_1 = 0.0;
                    double GBetaD_2 = 0.0;
                    double DBetaL_1 = 0.0;
                    double DBetaL_2 = 0.0;
                    double DBetaC = 0.0;
                    double DBetaG = 0.0;
                    double DBetaD = 0.0;
                    for (unsigned int f = 0; f < Dim; f++){
                        if (d == e){
                            GBetaA += tau_one(d,d) * fluid_fraction * viscosity * AGradN[j] * fluid_fraction_gradient[f] * rData.DN_DX(i,f);
                            AGBeta += tau_one(d,d) * fluid_fraction * viscosity * AGradN[i] * fluid_fraction_gradient[f] * rData.DN_DX(j,f);
                            GS += fluid_fraction * viscosity * rData.DN_DX(i,f) * rData.DN_DX(j,f);
                            AL += tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[j](f,f) * AGradN[i];
                            LA += tau_one(d,d) * std::pow(fluid_fraction,2) * viscosity * rData.DDN_DDX[i](f,f) * AGradN[j];
                            LL_diag_1 += tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](f,f);
                            LL_diag_2 += rData.DDN_DDX[j](f,f);
                            LGBeta_2 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DDN_DDX[i](f,f);
                            LGBeta_3 += rData.DN_DX(j,f) * fluid_fraction_gradient[f];
                            GBetaG_4 += tau_one(d,d) * std::pow(viscosity,2) * fluid_fraction_gradient[f] * rData.DN_DX(i,f);
                            GBetaG_5 += rData.DN_DX(j,f) * fluid_fraction_gradient[f];
                            GBetaL_4 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[f] * rData.DN_DX(i,f);
                            GBetaL_5 += rData.DDN_DDX[j](f,f);
                        }
                        LL_2 += tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](f,f) * rData.DDN_DDX[j](d,e);
                        LL_3 += tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](d,e) * rData.DDN_DDX[j](f,f);
                        LL_4 += tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](d,f) * rData.DDN_DDX[j](e,f);
                        LC_1 += 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](f,f) * rData.DDN_DDX[j](d,e);
                        LC_2 += 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](d,f) * rData.DDN_DDX[j](f,e);
                        CL_1 += 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](d,e) * rData.DDN_DDX[j](f,f);
                        CL_2 += 2.0 / 3.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](d,f) * rData.DDN_DDX[j](f,e);
                        CC += 4.0 / 9.0 * tau_one(d,d) * std::pow(fluid_fraction,2) * std::pow(viscosity,2) * rData.DDN_DDX[i](f,d) * rData.DDN_DDX[j](f,e);
                        LGBeta_1 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DDN_DDX[i](f,f) * rData.DN_DX(j,d) * fluid_fraction_gradient[e];
                        LGBeta_4 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DDN_DDX[i](d,f) * rData.DN_DX(j,f) * fluid_fraction_gradient[e];
                        LGBeta_5 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DDN_DDX[i](d,e) * fluid_fraction_gradient[f] * rData.DN_DX(j,f);
                        LDBeta_1 += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DN_DX(j,e) * fluid_fraction_gradient[d] * rData.DDN_DDX[i](f,f);
                        LDBeta_2 += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DN_DX(j,e) * fluid_fraction_gradient[f] * rData.DDN_DDX[i](d,f);
                        CGBeta_1 += 2.0 /3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DDN_DDX[i](d,f) * rData.DN_DX(j,f) * fluid_fraction_gradient[e];
                        CGBeta_2 += 2.0 /3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DDN_DDX[i](d,e) * rData.DN_DX(j,f) * fluid_fraction_gradient[f];
                        CDBeta += 4.0 / 9.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[f] * rData.DDN_DDX[i](d,f) * rData.DN_DX(j,e);
                        GBetaL_1 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(i,f) * rData.DDN_DDX[j](f,e);
                        GBetaL_2 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(i,e) * rData.DDN_DDX[j](f,f);
                        GBetaL_3 += tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[f] * rData.DN_DX(i,f) * rData.DDN_DDX[j](d,e);
                        GBetaC_1 += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(i,f) * rData.DDN_DDX[j](f,e);
                        GBetaC_2 += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * fluid_fraction_gradient[f] * rData.DN_DX(i,f) * rData.DDN_DDX[j](d,e);
                        DBetaL_1 += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.DDN_DDX[j](f,f);
                        DBetaL_2 += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DN_DX(i,d) * (fluid_fraction_gradient[f] * rData.DDN_DDX[j](f,e));
                        DBetaC += 4.0 / 9.0 * tau_one(d,d) * fluid_fraction * std::pow(viscosity,2) * rData.DN_DX(i,d) * fluid_fraction_gradient[f] * rData.DDN_DDX[j](e,f);
                        DBetaG += 4.0 / 3.0 * tau_one(d,d) * std::pow(viscosity,2) * rData.DN_DX(i,d) * fluid_fraction_gradient[f] * rData.DN_DX(j,f) * fluid_fraction_gradient[e];
                        GBetaG_1 += tau_one(d,d) * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(i,f) * rData.DN_DX(j,f) * fluid_fraction_gradient[e];
                        GBetaG_2 += tau_one(d,d) * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(i,e) * rData.DN_DX(j,f) * fluid_fraction_gradient[f];
                        GBetaG_3 += tau_one(d,d) * std::pow(viscosity,2) * fluid_fraction_gradient[f] * rData.DN_DX(i,f) * rData.DN_DX(j,d) * fluid_fraction_gradient[e];
                        GBetaD_1 += 2.0 / 3.0 * tau_one(d,d) * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(j,e) * fluid_fraction_gradient[f] * rData.DN_DX(i,f);
                        GBetaD_2 += 2.0 / 3.0 * tau_one(d,d) * std::pow(viscosity,2) * fluid_fraction_gradient[d] * rData.DN_DX(j,e) * fluid_fraction_gradient[f] * rData.DN_DX(i,f);
                        DBetaD += 4.0 / 9.0 * tau_one(d,d) * std::pow(viscosity, 2) * fluid_fraction_gradient[f] * fluid_fraction_gradient[f] * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    }
                    double LL = (LL_diag_1 * LL_diag_2) + LL_2 + LL_3 + LL_4;
                    double LGBeta = LGBeta_1 + (LGBeta_2 * LGBeta_3) + LGBeta_4 + LGBeta_5;
                    double LC = LC_1 + LC_2;
                    double CL = CL_1 + CL_2;
                    double LDBeta = LDBeta_1 + LDBeta_2;
                    double CGBeta = CGBeta_1 + CGBeta_2;
                    double GBetaL = GBetaL_1 + GBetaL_2 + GBetaL_3 + (GBetaL_4 * GBetaL_5);
                    double GBetaC = GBetaC_1 + GBetaC_2;
                    double DBetaL = DBetaL_1 + DBetaL_2;
                    double GBetaG = GBetaG_1 + GBetaG_2 + GBetaG_3 + (GBetaG_4 * GBetaG_5);
                    double GBetaD = GBetaD_1 + GBetaD_2;

                    LHS(row+d,col+e) += rData.Weight * (GS - DnuD + DD + DU + GU + GD + GBetaA - GBetaG + GBetaD - DBetaA + DBetaG - DBetaD - AGBeta + ADBeta + RSigma - AL + LA + AC - CA + LC + CL - CC - LL - LGBeta + LDBeta + CGBeta - CDBeta - GBetaL + GBetaC + DBetaL - DBetaC);

                }
                double GGBeta = GGBeta_1 + GGBeta_2;
                double GG = GG_1 + GG_2;

                LHS(row+Dim,col+d) += rData.Weight * (GA - GL + GC + QD - GGBeta + GDBeta + GAlphaD);

                LHS(row+d,col+Dim) += rData.Weight * (AG + LG - CG - P - GP + GG - DG);

            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * G;
        }

        // RHS terms
        double QAlphaF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            // v*BodyForce
            double VF = rData.N[i] * body_force[d];
            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            double AF = tau_one(d,d) * fluid_fraction * AGradN[i] * (body_force[d] - MomentumProj[d]);
            double LF = 0.0;
            double CF = 0.0;
            double GBetaF = 0.0;
            double DBetaF = 0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                LF += tau_one(d,d) * fluid_fraction * viscosity * rData.DDN_DDX[i](d,e) * (body_force[e] - MomentumProj[e]);
                LF += tau_one(d,d) * fluid_fraction * viscosity * rData.DDN_DDX[i](e,e) * (body_force[d] - MomentumProj[d]);
                CF += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * viscosity * rData.DDN_DDX[i](d,e) * (body_force[e] - MomentumProj[e]);
                GBetaF += tau_one(d,d) * viscosity * fluid_fraction_gradient[d] * rData.DN_DX(i,e) * (body_force[e] - MomentumProj[e]);
                GBetaF += tau_one(d,d) * viscosity * fluid_fraction_gradient[e] * rData.DN_DX(i,e) * (body_force[d] - MomentumProj[d]);
                DBetaF += 2.0 / 3.0 * tau_one(d,d) * viscosity * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * (body_force[e] - MomentumProj[e]);
            }
            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            QAlphaF += tau_one(d,d) * rData.DN_DX(i,d) * fluid_fraction * (body_force[d] - MomentumProj[d]);
            // OSS pressure subscale projection
            double VPhi = tau_two * rData.N[i] * fluid_fraction_gradient[d] * (mass_source - fluid_fraction_rate - MassProj);
            double DPhi = tau_two * rData.DN_DX(i,d) * fluid_fraction * (mass_source - fluid_fraction_rate - MassProj);
            rLocalRHS[row+d] += rData.Weight * (VF + AF + LF - CF + DPhi + GBetaF - DBetaF + VPhi);
        }
        double Q = rData.N[i] * (mass_source - fluid_fraction_rate);
        rLocalRHS[row+Dim] += rData.Weight * (QAlphaF + Q); // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Adding reactive terms to the stabilization
    if(!rData.UseOSS)
        this->AddReactionStabilization(rData,LHS,rLocalRHS);

    // Write (the linearized part of the) local contribution into residual form (A*dx = b - A*x)
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);
    noalias(rLocalRHS) -= prod(LHS, values);

    noalias(rLocalLHS) += LHS;
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        this->GetShapeSecondDerivatives(shape_function_second_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g],shape_function_second_derivatives[g]);

            this->AddMassLHS(data, rMassMatrix);
        }
    }
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                                                  VectorType& rRightHandSideVector,
                                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rDampMatrix.size1() != LocalSize )
        rDampMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        this->GetShapeSecondDerivatives(shape_function_second_derivatives);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            const auto& r_dndx = shape_derivatives[g];
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),r_dndx, shape_function_second_derivatives[g]);
            this->AddVelocitySystem(data, rDampMatrix, rRightHandSideVector);
        }
    }
}

template < class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateResistanceTensor(
    const TElementData& rData)
{
    BoundedMatrix<double,Dim,Dim>& rsigma = mViscousResistanceTensor[rData.IntegrationPointIndex];
    // typename TElementData::NodalVectorData variable_db;
    // const auto& rGeom = this->GetGeometry();
    // for (unsigned int i = 0; i < NumNodes; i++)
    //     for (unsigned int d = 0; d < Dim; d++)
    //         variable_db(i,d) = rGeom[i].FastGetSolutionStepValue(ACCELERATION)[d];
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);

    // const auto imposed_velocity = this->GetAtCoordinate(rData.Velocity, rData.N);
    // //const auto imposed_velocity = this->GetAtCoordinate(variable_db, rData.N);
    // const double porosity = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    // const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    // double kappa = (1 - porosity)/porosity;
    // //double kappa = 1.0;
    // double velocity_modulus = 0.0;

    // for (unsigned int d = 0; d < Dim; d++)
    //     velocity_modulus += imposed_velocity[d] * imposed_velocity[d];

    // double velocity_norm = std::sqrt(velocity_modulus);
    // for (unsigned int d = 0; d < Dim; d++)
        //rsigma(d,d) = (viscosity * 150.0 * std::pow(kappa,2) + 1.75 * kappa * velocity_norm);
    rsigma = permeability;
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::AddMassLHS(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;
            const double Mij = rData.Weight * density * fluid_fraction * rData.N[i] * rData.N[j];
            for (unsigned int d = 0; d < Dim; d++)
                rMassMatrix(row+d,col+d) += Mij;
        }
    }

    /* Note on OSS and full projection: Riccardo says that adding the terms provided by
     * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
     * could help reduce the non-linearity of the coupling between projection and u,p
     * However, leaving them on gives a lot of trouble whith the Bossak scheme:
     * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
     * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
     */
    if (!rData.UseOSS)
        this->AddMassStabilization(rData,rMassMatrix);
}

template<class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::MassProjTerm(
    const TElementData& rData,
    double &rMassRHS) const
{
        const auto velocities = rData.Velocity;
        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
        const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
        const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
        array_1d<double,Dim> fluid_fraction_gradient = ZeroVector(Dim);
        for (unsigned int i = 0; i < NumNodes; i++)
            for (unsigned int d = 0; d < Dim; d++)
                fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];

        // Compute this node's contribution to the residual (evaluated at integration point)
        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int d = 0; d < Dim; d++){
                rMassRHS -= (fluid_fraction_gradient[d] * rData.N[i] * velocities(i,d) + fluid_fraction * rData.DN_DX(i,d) * velocities(i,d));
            }
        }
        rMassRHS += mass_source - fluid_fraction_rate;
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateTau(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    BoundedMatrix<double,Dim,Dim> &TauOne,
    double &TauTwo) const
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;
    const int p = mInterpolationOrder;
    double inv_tau;
    double inv_tau_NS;
    //double spectral_radius = 0.0;
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction,rData.N);
    MatrixType sigma = ZeroMatrix(Dim+1, Dim+1);
    BoundedMatrix<double,Dim,Dim> I = IdentityMatrix(Dim, Dim);
    array_1d<double,Dim> fluid_fraction_gradient = ZeroVector(Dim);

    for (unsigned int i = 0; i < NumNodes; i++)
        for (unsigned int d = 0; d < Dim; d++)
            fluid_fraction_gradient[d] += rData.DN_DX(i,d) * rData.FluidFraction[i];


    double velocity_modulus = 0.0;
    double fluid_fraction_gradient_modulus = 0.0;

    for (unsigned int d = 0; d < Dim; d++){
        velocity_modulus += Velocity[d] * Velocity[d];
        fluid_fraction_gradient_modulus += std::pow(fluid_fraction_gradient[d],2);
        sigma(d,d) = mViscousResistanceTensor[rData.IntegrationPointIndex](d,d);
    }

    // This last term does not exist physically and it is included to do the spectral radius taking into account the inverse Gamma
    // whose size is (d+1,d+1)

    double velocity_norm = std::sqrt(velocity_modulus);
    double fluid_fraction_gradient_norm = std::sqrt(fluid_fraction_gradient_modulus);

    double c_alpha = fluid_fraction + h / c1 * fluid_fraction_gradient_norm;

    inv_tau_NS = c1 * viscosity / std::pow(h/(p*p),2.0) + density * (c2 * velocity_norm / (h/p) );
    double tau_one_NS = 1.0 / inv_tau_NS;

    //this->CalculateSpectralRadius(rData, spectral_radius, tau_one_NS, c1, sigma);

    inv_tau = inv_tau_NS * c_alpha;
    double tau_one = 1.0 / (inv_tau + sigma(0,0));

    TauOne = tau_one * I;
    TauTwo = std::pow(h/p,2.0) / (c1 * fluid_fraction * tau_one_NS);
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{
    // Get Shape function data

    DenseVector<DenseVector<Matrix>> ShapeFunctionSecondDerivatives;
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const long unsigned int NumGauss = GaussWeights.size();

    GeometryType& r_geometry = this->GetGeometry();
    this->GetShapeSecondDerivatives(ShapeFunctionSecondDerivatives);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    array_1d<double,NumNodes*Dim> momentum_rhs = ZeroVector(NumNodes*Dim);
    VectorType MassRHS = ZeroVector(NumNodes);
    VectorType NodalArea = ZeroVector(NumNodes);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        this->UpdateIntegrationPointData(data, g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g],ShapeFunctionSecondDerivatives[g]);

        array_1d<double, 3> MomentumRes = ZeroVector(3);
        double MassRes = 0.0;
        //array_1d<double,3> convective_velocity = this->GetAtCoordinate(data.Velocity,data.N) - this->GetAtCoordinate(data.MeshVelocity,data.N);
        array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(data);
        this->MomentumProjTerm(data, convective_velocity, MomentumRes);
        this->MassProjTerm(data,MassRes);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double W = data.Weight*data.N[i];
            unsigned int row = i*Dim;
            for (unsigned int d = 0; d < Dim; d++)
                momentum_rhs[row+d] += W*MomentumRes[d];
            MassRHS[i] += W*MassRes;
            NodalArea[i] += data.Weight*data.N[i];
        }
    }

    // Add carefully to nodal variables to avoid OpenMP race condition
    for (SizeType i = 0; i < NumNodes; ++i)
    {
        r_geometry[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double,3>& rMomValue = r_geometry[i].FastGetSolutionStepValue(ADVPROJ);
        unsigned int row = i*Dim;
        for (unsigned int d = 0; d < Dim; d++)
            rMomValue[d] += momentum_rhs[row+d];
        r_geometry[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
        r_geometry[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
        r_geometry[i].UnSetLock(); // Free the node for other threads
    }
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateSpectralRadius(
    const TElementData& rData,
    double& spectral_radius,
    double tau_one_NS,
    const double c1,
    MatrixType matrix) const
{
    const double h = rData.ElementSize;
    MatrixType inv_Gamma = IdentityMatrix(Dim+1, Dim+1);
    MatrixType eigen_val_mat, resulting_mat;
    MatrixType eigen_vect_mat = IdentityMatrix(Dim+1, Dim+1);

    double ji = h*h /(c1*std::pow(tau_one_NS,2));

    inv_Gamma(Dim,Dim) = 1/ji;
    resulting_mat = prod(inv_Gamma, eigen_vect_mat);

    this->GaussSeidelEigenSystem(matrix, resulting_mat, eigen_val_mat);

    for(unsigned int d = 0; d <= Dim; d++)
        for(unsigned int e = 0; e <= Dim; e++)
            spectral_radius = std::max(eigen_val_mat(d,e), spectral_radius);
}

template< class TElementData >
bool AlternativeQSVMSDEMCoupled<TElementData>::GaussSeidelEigenSystem(
        MatrixType& rA,
        MatrixType& rEigenVectorsMatrix,
        MatrixType& rEigenValuesMatrix,
        const double Tolerance,
        const SizeType MaxIterations
        ) const
    {
        bool is_converged = false;

        const SizeType size = rA.size1();

#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        KRATOS_WARNING_IF("EigenSystem", rEigenVectorsMatrix.size1() != size || rEigenVectorsMatrix.size2() != size) << "EigenSystem has detected an incorrect size of your EigenVectorsMatrix matrix. Please resize before compute" << std::endl;
        KRATOS_WARNING_IF("EigenSystem", rEigenValuesMatrix.size1() != size || rEigenValuesMatrix.size2() != size) << "EigenSystem has detected an incorrect size of your EigenValuesMatrix matrix. Please resize before compute" << std::endl;
#else
        if (rEigenVectorsMatrix.size1() != size || rEigenVectorsMatrix.size2() != size)
            rEigenVectorsMatrix.resize(size, size, false);
        if (rEigenValuesMatrix.size1() != size || rEigenValuesMatrix.size2() != size)
            rEigenValuesMatrix.resize(size, size, false);
#endif // KRATOS_USE_AMATRIX

        const MatrixType identity_matrix = IdentityMatrix(size);
        noalias(rEigenValuesMatrix) = rA;

        // Auxiliar values
        MatrixType aux_A, aux_V_matrix, rotation_matrix;
        double a, u, c, s, gamma, teta;
        IndexType index1, index2;

        aux_A.resize(size,size,false);
        aux_V_matrix.resize(size,size,false);
        rotation_matrix.resize(size,size,false);

        for(IndexType iterations = 0; iterations < MaxIterations; ++iterations) {
            is_converged = true;

            a = 0.0;
            index1 = 0;
            index2 = 1;

            for(IndexType i = 0; i < size; ++i) {
                for(IndexType j = (i + 1); j < size; ++j) {
                    if((std::abs(rEigenValuesMatrix(i, j)) > a ) && (std::abs(rEigenValuesMatrix(i, j)) > Tolerance)) {
                        a = std::abs(rEigenValuesMatrix(i,j));
                        index1 = i;
                        index2 = j;
                        is_converged = false;
                    }
                }
            }

            if(is_converged) {
                break;
            }

            // Calculation of Rotation angle
            gamma = (rEigenValuesMatrix(index2, index2)-rEigenValuesMatrix(index1, index1)) / (2 * rEigenValuesMatrix(index1, index2));
            u = 1.0;

            if(std::abs(gamma) > Tolerance && std::abs(gamma)< (1.0/Tolerance)) {
                u = gamma / std::abs(gamma) * 1.0 / (std::abs(gamma) + std::sqrt(1.0 + gamma * gamma));
            } else {
                if  (std::abs(gamma) >= (1.0/Tolerance)) {
                    u = 0.5 / gamma;
                }
            }

            c = 1.0 / (std::sqrt(1.0 + u * u));
            s = c * u;
            teta = s / (1.0 + c);

            // Rotation of the Matrix
            noalias(aux_A) = rEigenValuesMatrix;
            aux_A(index2, index2) = rEigenValuesMatrix(index2,index2) + u * rEigenValuesMatrix(index1, index2);
            aux_A(index1, index1) = rEigenValuesMatrix(index1,index1) - u * rEigenValuesMatrix(index1, index2);
            aux_A(index1, index2) = 0.0;
            aux_A(index2, index1) = 0.0;

            for(IndexType i = 0; i < size; ++i) {
                if((i!= index1) && (i!= index2)) {
                    aux_A(index2, i) = rEigenValuesMatrix(index2, i) + s * (rEigenValuesMatrix(index1, i)- teta * rEigenValuesMatrix(index2, i));
                    aux_A(i, index2) = rEigenValuesMatrix(index2, i) + s * (rEigenValuesMatrix(index1, i)- teta * rEigenValuesMatrix(index2, i));
                    aux_A(index1, i) = rEigenValuesMatrix(index1, i) - s * (rEigenValuesMatrix(index2, i) + teta * rEigenValuesMatrix(index1, i));
                    aux_A(i, index1) = rEigenValuesMatrix(index1, i) - s * (rEigenValuesMatrix(index2, i) + teta * rEigenValuesMatrix(index1, i));
                }
            }

            noalias(rEigenValuesMatrix) = aux_A;

            // Calculation of the eigeneigen vector matrix V
            noalias(rotation_matrix) = identity_matrix;
            rotation_matrix(index2, index1) = -s;
            rotation_matrix(index1, index2) =  s;
            rotation_matrix(index1, index1) =  c;
            rotation_matrix(index2, index2) =  c;

            noalias(aux_V_matrix) = ZeroMatrix(size, size);

            for(IndexType i = 0; i < size; ++i) {
                for(IndexType j = 0; j < size; ++j) {
                    for(IndexType k = 0; k < size; ++k) {
                        aux_V_matrix(i, j) += rEigenVectorsMatrix(i, k) * rotation_matrix(k, j);
                    }
                }
            }
            noalias(rEigenVectorsMatrix) = aux_V_matrix;
        }

        KRATOS_WARNING_IF("MathUtils::EigenSystem", !is_converged) << "Spectral decomposition not converged " << std::endl;

        return is_converged;
    }

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    array_1d<double,3> &rVelocitySubscale) const
{
    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;

    //array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    array_1d<double,3> Residual = ZeroVector(3);

    if (!rData.UseOSS)
        this->AlgebraicMomentumResidual(rData,convective_velocity,Residual);
    else
        this->OrthogonalMomentumResidual(rData,convective_velocity,Residual);

    for (unsigned int d = 0; d < Dim; ++d)
        rVelocitySubscale[d] = tau_one(d,d) * Residual[d];
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::SubscalePressure(
        const TElementData& rData,
        double &rPressureSubscale) const
{
    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    //const auto& rGeom = this->GetGeometry();
    //typename TElementData::NodalVectorData variable_db;
    // for (unsigned int i = 0; i < NumNodes; i++)
    //     for (unsigned int d = 0; d < Dim; d++)
    //         variable_db(i,d) = rGeom[i].FastGetSolutionStepValue(ACCELERATION)[d];
    // const array_1d<double, 3> convective_velocity =
    //     this->GetAtCoordinate(variable_db, rData.N) -
    //     this->GetAtCoordinate(rData.MeshVelocity, rData.N);
    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);
    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    double Residual = 0.0;

    if (!rData.UseOSS)
        this->AlgebraicMassResidual(rData,Residual);
    else
        this->OrthogonalMassResidual(rData,Residual);

    rPressureSubscale = tau_two*Residual;
}

template< class TElementData >
array_1d<double,3> AlternativeQSVMSDEMCoupled<TElementData>::FullConvectiveVelocity(
    const TElementData& rData) const
{
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    //Adding subscale term componentwise because return type is of size 3, but subscale is of size Dim
    const array_1d<double,Dim>& r_predicted_subscale = mPredictedSubscaleVelocity[rData.IntegrationPointIndex];

    for (unsigned int d = 0; d < Dim; d++) {
        convective_velocity[d] += r_predicted_subscale[d];
    }

    return convective_velocity;
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::UpdateSubscaleVelocity(
    const TElementData& rData)
{
    array_1d<double,Dim> predicted_subscale_velocity;

    array_1d<double,Dim> previous_velocity = mPreviousVelocity[rData.IntegrationPointIndex];

    //for (size_t i = 0; i < NumNodes; i++) {
    array_1d<double,Dim> subscale_velocity_on_previous_iteration = mPredictedSubscaleVelocity[rData.IntegrationPointIndex];

    array_1d<double,3> v_d = ZeroVector(3);

    for (unsigned int d = 0; d < Dim; d++)
        v_d[d] = previous_velocity[d] + subscale_velocity_on_previous_iteration[d];

    // Part of the residual that does not depend on the subscale
    array_1d<double,3> static_residual = ZeroVector(3);

    if (!rData.UseOSS)
        this->AlgebraicMomentumResidual(rData,v_d,static_residual);
    else
        this->OrthogonalMomentumResidual(rData,v_d,static_residual);


    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;

    this->CalculateTau(rData,v_d,tau_one,tau_two);

    for (unsigned int d = 0; d < Dim; d++)
        predicted_subscale_velocity[d] = tau_one(d,d) * static_residual[d];

    noalias(mPredictedSubscaleVelocity[rData.IntegrationPointIndex]) = predicted_subscale_velocity;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::save(Serializer& rSerializer) const
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<2,3> >;
template class AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<3,4> >;

template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<2,4> >;
template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<2,6> >;
template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<2,9> >;

template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;
template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<3,27> >;
}

