//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Uxue Chasco
//

#include "vectorial_convection_fractional_element.h"

namespace Kratos
{

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    VectorialConvectionFractionalElement<TElementData>::VectorialConvectionFractionalElement(IndexType NewId) : Element(NewId)
    {
    }

    template <class TElementData>
    VectorialConvectionFractionalElement<TElementData>::VectorialConvectionFractionalElement(IndexType NewId, const NodesArrayType &ThisNodes) : Element(NewId, ThisNodes)
    {
    }

    template <class TElementData>
    VectorialConvectionFractionalElement<TElementData>::VectorialConvectionFractionalElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
    }

    template <class TElementData>
    VectorialConvectionFractionalElement<TElementData>::VectorialConvectionFractionalElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties)
    {
    }

    template <class TElementData>
    VectorialConvectionFractionalElement<TElementData>::~VectorialConvectionFractionalElement()
    {
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    Element::Pointer VectorialConvectionFractionalElement<TElementData>::Create(
        IndexType NewId,
        NodesArrayType const &ThisNodes,
        Properties::Pointer pProperties) const
    {
        return Kratos::make_intrusive<VectorialConvectionFractionalElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /***********************************************************************************/

    template <class TElementData>
    Element::Pointer VectorialConvectionFractionalElement<TElementData>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const
    {
        return Kratos::make_intrusive<VectorialConvectionFractionalElement>(NewId, pGeom, pProperties);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::CalculateRightHandSide(VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented", "");
    }

    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();
        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        {
            UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::UpdateIntegrationPointData(
        TElementData &rData,
        unsigned int IntegrationPointIndex,
        double Weight,
        const typename TElementData::MatrixRowType &rN,
        const typename TElementData::ShapeDerivativesType &rDN_DX) const
    {
        rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::AddTimeIntegratedSystem(
        TElementData &rData,
        MatrixType &rLHS,
        VectorType &rRHS)
    {
        this->ComputeGaussPointLHSContribution(rData, rLHS);
        this->ComputeGaussPointRHSContribution(rData, rRHS);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
                                                                                   Matrix &rNContainer,
                                                                                   ShapeFunctionDerivativesArrayType &rDN_DX) const
    {
        const GeometryData::IntegrationMethod integration_method = this->GetIntegrationMethod();
        const GeometryType &r_geometry = this->GetGeometry();
        const unsigned int number_of_gauss_points = r_geometry.IntegrationPointsNumber(integration_method);

        Vector DetJ;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, integration_method);

        if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != NumNodes)
        {
            rNContainer.resize(number_of_gauss_points, NumNodes, false);
        }
        rNContainer = r_geometry.ShapeFunctionsValues(integration_method);

        const GeometryType::IntegrationPointsArrayType &IntegrationPoints = r_geometry.IntegrationPoints(integration_method);

        if (rGaussWeights.size() != number_of_gauss_points)
        {
            rGaussWeights.resize(number_of_gauss_points, false);
        }

        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<2, 3>>::ComputeGaussPointLHSContribution(VectorialConvectionFractionalElementData<2, 3> &rData,
                                                                                                                                MatrixType &rLHS)
    {

        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double dyn_tau = rData.DynamicTau;
        const auto &vfrac = rData.Velocity_Fractional;
        const auto &vmesh = rData.MeshVelocity;
        const auto &vconv = vfrac-vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c2 = 2.0;

        auto &lhs = rData.lhs;
        //substitute_lhs_2D
        //  Add intermediate results to local system.
        noalias(rLHS) += lhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<3, 4>>::ComputeGaussPointLHSContribution(VectorialConvectionFractionalElementData<3, 4> &rData,
                                                                                                                                MatrixType &rLHS)
    {


        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;

        const double dyn_tau = rData.DynamicTau;
        const auto &vfrac = rData.Velocity_Fractional;
        const auto &vmesh = rData.MeshVelocity;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c2 = 2.0;

        auto &lhs = rData.lhs;

        //substitute_lhs_3D
        //  Add intermediate results to local system.
        noalias(rLHS) += lhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<2, 3>>::ComputeGaussPointRHSContribution(VectorialConvectionFractionalElementData<2, 3> &rData,
                                                                                                                                VectorType &rRHS)
    {

        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double bdf1 = rData.bdf1;
        const double bdf2 = rData.bdf2;
        const double dyn_tau = rData.DynamicTau;
        const auto &vn = rData.Velocity_OldStep1;
        const auto &vnn = rData.Velocity_OldStep2;
        // const auto &vnnn = rData.Velocity_OldStep3; #an bdf2
        const auto &vmesh = rData.MeshVelocity;
        const auto &vfrac = rData.Velocity_Fractional;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c2 = 2.0;

        auto &rhs = rData.rhs;

        //substitute_rhs_2D
        //Add intermediate results to local system.
        noalias(rRHS) += rhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<3, 4>>::ComputeGaussPointRHSContribution(VectorialConvectionFractionalElementData<3, 4> &rData,
                                                                                                                                VectorType &rRHS)
    {


        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double bdf1 = rData.bdf1;
        const double bdf2 = rData.bdf2;
        const double dyn_tau = rData.DynamicTau;
        const auto &vn = rData.Velocity_OldStep1;
        const auto &vnn = rData.Velocity_OldStep2;
        // const auto &vnnn = rData.Velocity_OldStep3; #an_bdf2
        const auto &vmesh = rData.MeshVelocity;
        const auto &vfrac = rData.Velocity_Fractional;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c2 = 2.0;

        auto &rhs = rData.rhs;

        //substitute_rhs_3D

        //  Add intermediate results to local system.
        noalias(rRHS) += rhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const GeometryType &r_geometry = this->GetGeometry();

        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACTIONAL_VELOCITY_X);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(FRACTIONAL_VELOCITY_X, xpos).EquationId();
            rResult[LocalIndex++] = r_geometry[i].GetDof(FRACTIONAL_VELOCITY_Y, xpos + 1).EquationId();
            if (Dim == 3)
                rResult[LocalIndex++] = r_geometry[i].GetDof(FRACTIONAL_VELOCITY_Z, xpos + 2).EquationId();
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    void VectorialConvectionFractionalElement<TElementData>::GetDofList(
        DofsVectorType &rElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const GeometryType &r_geometry = this->GetGeometry();

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACTIONAL_VELOCITY_X);

        unsigned int LocalIndex = 0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(FRACTIONAL_VELOCITY_X, xpos);
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(FRACTIONAL_VELOCITY_Y, xpos + 1);
            if (Dim == 3)
                rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(FRACTIONAL_VELOCITY_Z, xpos + 2);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /***********************************************************************************/
    /***********************************************************************************/
    template class VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<2, 3>>;
    template class VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<3, 4>>;

    /***********************************************************************************/
    /***********************************************************************************/
}
