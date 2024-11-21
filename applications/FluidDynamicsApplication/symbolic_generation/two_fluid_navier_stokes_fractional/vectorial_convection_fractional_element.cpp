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

        const double mu = rData.EffectiveViscosity;
        const double h = rData.ElementSize;

        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;

        const double dyn_tau = rData.DynamicTau;

        const auto vfrac = rData.Velocity_Fractional;
        const auto vn = rData.Velocity_OldStep1;


        const auto vmesh = rData.MeshVelocity;
        const auto vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c1 = 4.0;
        constexpr double stab_c2 = 2.0;

        auto &lhs = rData.lhs;
        const double clhs0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs2 = DN(0,0)*clhs0 + DN(0,1)*clhs1;
const double clhs3 = N[0]*bdf0;
const double clhs4 = clhs2 + clhs3;
const double clhs5 = 1.0/(stab_c2*sqrt(pow(clhs0, 2) + pow(clhs1, 2))/h + dyn_tau/dt);
const double clhs6 = clhs5*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double clhs7 = N[0]*clhs6;
const double clhs8 = clhs2*clhs5;
const double clhs9 = pow(N[0], 2)*bdf0 + N[0]*clhs2 + clhs4*clhs7 + clhs4*clhs8;
const double clhs10 = N[1]*clhs3;
const double clhs11 = DN(1,0)*clhs0 + DN(1,1)*clhs1;
const double clhs12 = N[1]*bdf0;
const double clhs13 = clhs11 + clhs12;
const double clhs14 = N[0]*clhs11 + clhs10 + clhs13*clhs7 + clhs13*clhs8;
const double clhs15 = N[2]*clhs3;
const double clhs16 = DN(2,0)*clhs0 + DN(2,1)*clhs1;
const double clhs17 = N[2]*bdf0 + clhs16;
const double clhs18 = N[0]*clhs16 + clhs15 + clhs17*clhs7 + clhs17*clhs8;
const double clhs19 = N[1]*clhs6;
const double clhs20 = clhs11*clhs5;
const double clhs21 = N[1]*clhs2 + clhs10 + clhs19*clhs4 + clhs20*clhs4;
const double clhs22 = pow(N[1], 2)*bdf0 + N[1]*clhs11 + clhs13*clhs19 + clhs13*clhs20;
const double clhs23 = N[2]*clhs12;
const double clhs24 = N[1]*clhs16 + clhs17*clhs19 + clhs17*clhs20 + clhs23;
const double clhs25 = N[2]*clhs6;
const double clhs26 = clhs16*clhs5;
const double clhs27 = N[2]*clhs2 + clhs15 + clhs25*clhs4 + clhs26*clhs4;
const double clhs28 = N[2]*clhs11 + clhs13*clhs25 + clhs13*clhs26 + clhs23;
const double clhs29 = pow(N[2], 2)*bdf0 + N[2]*clhs16 + clhs17*clhs25 + clhs17*clhs26;
lhs(0,0)=clhs9;
lhs(0,1)=0;
lhs(0,2)=clhs14;
lhs(0,3)=0;
lhs(0,4)=clhs18;
lhs(0,5)=0;
lhs(1,0)=0;
lhs(1,1)=clhs9;
lhs(1,2)=0;
lhs(1,3)=clhs14;
lhs(1,4)=0;
lhs(1,5)=clhs18;
lhs(2,0)=clhs21;
lhs(2,1)=0;
lhs(2,2)=clhs22;
lhs(2,3)=0;
lhs(2,4)=clhs24;
lhs(2,5)=0;
lhs(3,0)=0;
lhs(3,1)=clhs21;
lhs(3,2)=0;
lhs(3,3)=clhs22;
lhs(3,4)=0;
lhs(3,5)=clhs24;
lhs(4,0)=clhs27;
lhs(4,1)=0;
lhs(4,2)=clhs28;
lhs(4,3)=0;
lhs(4,4)=clhs29;
lhs(4,5)=0;
lhs(5,0)=0;
lhs(5,1)=clhs27;
lhs(5,2)=0;
lhs(5,3)=clhs28;
lhs(5,4)=0;
lhs(5,5)=clhs29;

        //  Add intermediate results to local system.
        noalias(rLHS) += lhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<3, 4>>::ComputeGaussPointLHSContribution(VectorialConvectionFractionalElementData<3, 4> &rData,
                                                                                                                                MatrixType &rLHS)
    {

        const double mu = rData.EffectiveViscosity;
        const double h = rData.ElementSize;

        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;

        const double dyn_tau = rData.DynamicTau;
        const auto vn = rData.Velocity_OldStep1;
        const auto vfrac = rData.Velocity_Fractional;
        const auto vmesh = rData.MeshVelocity;
        const auto vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c1 = 4.0;
        constexpr double stab_c2 = 2.0;

        auto &lhs = rData.lhs;

        const double clhs0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs3 = DN(0,0)*clhs0 + DN(0,1)*clhs1 + DN(0,2)*clhs2;
const double clhs4 = N[0]*bdf0;
const double clhs5 = clhs3 + clhs4;
const double clhs6 = 1.0/(stab_c2*sqrt(pow(clhs0, 2) + pow(clhs1, 2) + pow(clhs2, 2))/h + dyn_tau/dt);
const double clhs7 = clhs6*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double clhs8 = N[0]*clhs7;
const double clhs9 = clhs3*clhs6;
const double clhs10 = pow(N[0], 2)*bdf0 + N[0]*clhs3 + clhs5*clhs8 + clhs5*clhs9;
const double clhs11 = N[1]*clhs4;
const double clhs12 = DN(1,0)*clhs0 + DN(1,1)*clhs1 + DN(1,2)*clhs2;
const double clhs13 = N[1]*bdf0;
const double clhs14 = clhs12 + clhs13;
const double clhs15 = N[0]*clhs12 + clhs11 + clhs14*clhs8 + clhs14*clhs9;
const double clhs16 = N[2]*clhs4;
const double clhs17 = DN(2,0)*clhs0 + DN(2,1)*clhs1 + DN(2,2)*clhs2;
const double clhs18 = N[2]*bdf0;
const double clhs19 = clhs17 + clhs18;
const double clhs20 = N[0]*clhs17 + clhs16 + clhs19*clhs8 + clhs19*clhs9;
const double clhs21 = N[3]*clhs4;
const double clhs22 = DN(3,0)*clhs0 + DN(3,1)*clhs1 + DN(3,2)*clhs2;
const double clhs23 = N[3]*bdf0 + clhs22;
const double clhs24 = N[0]*clhs22 + clhs21 + clhs23*clhs8 + clhs23*clhs9;
const double clhs25 = N[1]*clhs7;
const double clhs26 = clhs12*clhs6;
const double clhs27 = N[1]*clhs3 + clhs11 + clhs25*clhs5 + clhs26*clhs5;
const double clhs28 = pow(N[1], 2)*bdf0 + N[1]*clhs12 + clhs14*clhs25 + clhs14*clhs26;
const double clhs29 = N[2]*clhs13;
const double clhs30 = N[1]*clhs17 + clhs19*clhs25 + clhs19*clhs26 + clhs29;
const double clhs31 = N[3]*clhs13;
const double clhs32 = N[1]*clhs22 + clhs23*clhs25 + clhs23*clhs26 + clhs31;
const double clhs33 = N[2]*clhs7;
const double clhs34 = clhs17*clhs6;
const double clhs35 = N[2]*clhs3 + clhs16 + clhs33*clhs5 + clhs34*clhs5;
const double clhs36 = N[2]*clhs12 + clhs14*clhs33 + clhs14*clhs34 + clhs29;
const double clhs37 = pow(N[2], 2)*bdf0 + N[2]*clhs17 + clhs19*clhs33 + clhs19*clhs34;
const double clhs38 = N[3]*clhs18;
const double clhs39 = N[2]*clhs22 + clhs23*clhs33 + clhs23*clhs34 + clhs38;
const double clhs40 = N[3]*clhs7;
const double clhs41 = clhs22*clhs6;
const double clhs42 = N[3]*clhs3 + clhs21 + clhs40*clhs5 + clhs41*clhs5;
const double clhs43 = N[3]*clhs12 + clhs14*clhs40 + clhs14*clhs41 + clhs31;
const double clhs44 = N[3]*clhs17 + clhs19*clhs40 + clhs19*clhs41 + clhs38;
const double clhs45 = pow(N[3], 2)*bdf0 + N[3]*clhs22 + clhs23*clhs40 + clhs23*clhs41;
lhs(0,0)=clhs10;
lhs(0,1)=0;
lhs(0,2)=0;
lhs(0,3)=clhs15;
lhs(0,4)=0;
lhs(0,5)=0;
lhs(0,6)=clhs20;
lhs(0,7)=0;
lhs(0,8)=0;
lhs(0,9)=clhs24;
lhs(0,10)=0;
lhs(0,11)=0;
lhs(1,0)=0;
lhs(1,1)=clhs10;
lhs(1,2)=0;
lhs(1,3)=0;
lhs(1,4)=clhs15;
lhs(1,5)=0;
lhs(1,6)=0;
lhs(1,7)=clhs20;
lhs(1,8)=0;
lhs(1,9)=0;
lhs(1,10)=clhs24;
lhs(1,11)=0;
lhs(2,0)=0;
lhs(2,1)=0;
lhs(2,2)=clhs10;
lhs(2,3)=0;
lhs(2,4)=0;
lhs(2,5)=clhs15;
lhs(2,6)=0;
lhs(2,7)=0;
lhs(2,8)=clhs20;
lhs(2,9)=0;
lhs(2,10)=0;
lhs(2,11)=clhs24;
lhs(3,0)=clhs27;
lhs(3,1)=0;
lhs(3,2)=0;
lhs(3,3)=clhs28;
lhs(3,4)=0;
lhs(3,5)=0;
lhs(3,6)=clhs30;
lhs(3,7)=0;
lhs(3,8)=0;
lhs(3,9)=clhs32;
lhs(3,10)=0;
lhs(3,11)=0;
lhs(4,0)=0;
lhs(4,1)=clhs27;
lhs(4,2)=0;
lhs(4,3)=0;
lhs(4,4)=clhs28;
lhs(4,5)=0;
lhs(4,6)=0;
lhs(4,7)=clhs30;
lhs(4,8)=0;
lhs(4,9)=0;
lhs(4,10)=clhs32;
lhs(4,11)=0;
lhs(5,0)=0;
lhs(5,1)=0;
lhs(5,2)=clhs27;
lhs(5,3)=0;
lhs(5,4)=0;
lhs(5,5)=clhs28;
lhs(5,6)=0;
lhs(5,7)=0;
lhs(5,8)=clhs30;
lhs(5,9)=0;
lhs(5,10)=0;
lhs(5,11)=clhs32;
lhs(6,0)=clhs35;
lhs(6,1)=0;
lhs(6,2)=0;
lhs(6,3)=clhs36;
lhs(6,4)=0;
lhs(6,5)=0;
lhs(6,6)=clhs37;
lhs(6,7)=0;
lhs(6,8)=0;
lhs(6,9)=clhs39;
lhs(6,10)=0;
lhs(6,11)=0;
lhs(7,0)=0;
lhs(7,1)=clhs35;
lhs(7,2)=0;
lhs(7,3)=0;
lhs(7,4)=clhs36;
lhs(7,5)=0;
lhs(7,6)=0;
lhs(7,7)=clhs37;
lhs(7,8)=0;
lhs(7,9)=0;
lhs(7,10)=clhs39;
lhs(7,11)=0;
lhs(8,0)=0;
lhs(8,1)=0;
lhs(8,2)=clhs35;
lhs(8,3)=0;
lhs(8,4)=0;
lhs(8,5)=clhs36;
lhs(8,6)=0;
lhs(8,7)=0;
lhs(8,8)=clhs37;
lhs(8,9)=0;
lhs(8,10)=0;
lhs(8,11)=clhs39;
lhs(9,0)=clhs42;
lhs(9,1)=0;
lhs(9,2)=0;
lhs(9,3)=clhs43;
lhs(9,4)=0;
lhs(9,5)=0;
lhs(9,6)=clhs44;
lhs(9,7)=0;
lhs(9,8)=0;
lhs(9,9)=clhs45;
lhs(9,10)=0;
lhs(9,11)=0;
lhs(10,0)=0;
lhs(10,1)=clhs42;
lhs(10,2)=0;
lhs(10,3)=0;
lhs(10,4)=clhs43;
lhs(10,5)=0;
lhs(10,6)=0;
lhs(10,7)=clhs44;
lhs(10,8)=0;
lhs(10,9)=0;
lhs(10,10)=clhs45;
lhs(10,11)=0;
lhs(11,0)=0;
lhs(11,1)=0;
lhs(11,2)=clhs42;
lhs(11,3)=0;
lhs(11,4)=0;
lhs(11,5)=clhs43;
lhs(11,6)=0;
lhs(11,7)=0;
lhs(11,8)=clhs44;
lhs(11,9)=0;
lhs(11,10)=0;
lhs(11,11)=clhs45;

        //  Add intermediate results to local system.
        noalias(rLHS) += lhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<2, 3>>::ComputeGaussPointRHSContribution(VectorialConvectionFractionalElementData<2, 3> &rData,
                                                                                                                                VectorType &rRHS)
    {

        const double mu = rData.EffectiveViscosity;

        const double h = rData.ElementSize;

        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double bdf1 = rData.bdf1;
        const double bdf2 = rData.bdf2;

        const double dyn_tau = rData.DynamicTau;

        const auto &vn = rData.Velocity_OldStep1;

        const auto &vnnn = rData.Velocity_OldStep3;
        const auto &vnn = rData.Velocity_OldStep2;
        const auto &vmesh = rData.MeshVelocity;
        const auto vfrac = rData.Velocity_Fractional;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c1 = 4.0;
        constexpr double stab_c2 = 2.0;

        auto &rhs = rData.rhs;

        const double crhs0 = N[0]*(bdf0*vn(0,0) + bdf1*vnn(0,0) + bdf2*vnnn(0,0));
const double crhs1 = N[1]*(bdf0*vn(1,0) + bdf1*vnn(1,0) + bdf2*vnnn(1,0));
const double crhs2 = N[2]*(bdf0*vn(2,0) + bdf1*vnn(2,0) + bdf2*vnnn(2,0));
const double crhs3 = crhs0 + crhs1 + crhs2;
const double crhs4 = N[0]*(bdf0*vfrac(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*vfrac(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*vfrac(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
const double crhs5 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crhs6 = crhs5*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0));
const double crhs7 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crhs8 = crhs7*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0));
const double crhs9 = crhs6 + crhs8;
const double crhs10 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs11 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs12 = crhs10*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0)) + crhs11*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0));
const double crhs13 = -crhs0 - crhs1 + crhs12 - crhs2 + crhs4 - crhs6 - crhs8;
const double crhs14 = 1.0/(stab_c2*sqrt(pow(crhs10, 2) + pow(crhs11, 2))/h + dyn_tau/dt);
const double crhs15 = crhs14*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs16 = N[0]*crhs15;
const double crhs17 = crhs14*(DN(0,0)*crhs10 + DN(0,1)*crhs11);
const double crhs18 = N[0]*(bdf0*vn(0,1) + bdf1*vnn(0,1) + bdf2*vnnn(0,1));
const double crhs19 = N[1]*(bdf0*vn(1,1) + bdf1*vnn(1,1) + bdf2*vnnn(1,1));
const double crhs20 = N[2]*(bdf0*vn(2,1) + bdf1*vnn(2,1) + bdf2*vnnn(2,1));
const double crhs21 = crhs18 + crhs19 + crhs20;
const double crhs22 = N[0]*(bdf0*vfrac(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*vfrac(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*vfrac(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
const double crhs23 = crhs5*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1));
const double crhs24 = crhs7*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1));
const double crhs25 = crhs23 + crhs24;
const double crhs26 = crhs10*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1)) + crhs11*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1));
const double crhs27 = -crhs18 - crhs19 - crhs20 + crhs22 - crhs23 - crhs24 + crhs26;
const double crhs28 = N[1]*crhs15;
const double crhs29 = crhs14*(DN(1,0)*crhs10 + DN(1,1)*crhs11);
const double crhs30 = N[2]*crhs15;
const double crhs31 = crhs14*(DN(2,0)*crhs10 + DN(2,1)*crhs11);
rhs[0]=-N[0]*crhs12 + N[0]*crhs3 - N[0]*crhs4 + N[0]*crhs9 - crhs13*crhs16 - crhs13*crhs17;
rhs[1]=N[0]*crhs21 - N[0]*crhs22 + N[0]*crhs25 - N[0]*crhs26 - crhs16*crhs27 - crhs17*crhs27;
rhs[2]=-N[1]*crhs12 + N[1]*crhs3 - N[1]*crhs4 + N[1]*crhs9 - crhs13*crhs28 - crhs13*crhs29;
rhs[3]=N[1]*crhs21 - N[1]*crhs22 + N[1]*crhs25 - N[1]*crhs26 - crhs27*crhs28 - crhs27*crhs29;
rhs[4]=-N[2]*crhs12 + N[2]*crhs3 - N[2]*crhs4 + N[2]*crhs9 - crhs13*crhs30 - crhs13*crhs31;
rhs[5]=N[2]*crhs21 - N[2]*crhs22 + N[2]*crhs25 - N[2]*crhs26 - crhs27*crhs30 - crhs27*crhs31;

        //Add intermediate results to local system.
        noalias(rRHS) += rhs * rData.Weight;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void VectorialConvectionFractionalElement<VectorialConvectionFractionalElementData<3, 4>>::ComputeGaussPointRHSContribution(VectorialConvectionFractionalElementData<3, 4> &rData,
                                                                                                                                VectorType &rRHS)
    {

        const double mu = rData.EffectiveViscosity;

        const double h = rData.ElementSize;

        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double bdf1 = rData.bdf1;
        const double bdf2 = rData.bdf2;

        const double dyn_tau = rData.DynamicTau;

        const auto &vn = rData.Velocity_OldStep1;
        const auto &vnn = rData.Velocity_OldStep2;
        const auto &vnnn = rData.Velocity_OldStep3;

        const auto &vmesh = rData.MeshVelocity;
        const auto vfrac = rData.Velocity_Fractional;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c1 = 4.0;
        constexpr double stab_c2 = 2.0;

        auto &rhs = rData.rhs;

        const double crhs0 = N[0]*(bdf0*vn(0,0) + bdf1*vnn(0,0) + bdf2*vnnn(0,0));
const double crhs1 = N[1]*(bdf0*vn(1,0) + bdf1*vnn(1,0) + bdf2*vnnn(1,0));
const double crhs2 = N[2]*(bdf0*vn(2,0) + bdf1*vnn(2,0) + bdf2*vnnn(2,0));
const double crhs3 = N[3]*(bdf0*vn(3,0) + bdf1*vnn(3,0) + bdf2*vnnn(3,0));
const double crhs4 = crhs0 + crhs1 + crhs2 + crhs3;
const double crhs5 = N[0]*(bdf0*vfrac(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*vfrac(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*vfrac(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*vfrac(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
const double crhs6 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crhs7 = crhs6*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0));
const double crhs8 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crhs9 = crhs8*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0));
const double crhs10 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crhs11 = crhs10*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0));
const double crhs12 = crhs11 + crhs7 + crhs9;
const double crhs13 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs14 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs15 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs16 = crhs13*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0)) + crhs14*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0)) + crhs15*(DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0));
const double crhs17 = -crhs0 - crhs1 - crhs11 + crhs16 - crhs2 - crhs3 + crhs5 - crhs7 - crhs9;
const double crhs18 = 1.0/(stab_c2*sqrt(pow(crhs13, 2) + pow(crhs14, 2) + pow(crhs15, 2))/h + dyn_tau/dt);
const double crhs19 = crhs18*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs20 = N[0]*crhs19;
const double crhs21 = crhs18*(DN(0,0)*crhs13 + DN(0,1)*crhs14 + DN(0,2)*crhs15);
const double crhs22 = N[0]*(bdf0*vn(0,1) + bdf1*vnn(0,1) + bdf2*vnnn(0,1));
const double crhs23 = N[1]*(bdf0*vn(1,1) + bdf1*vnn(1,1) + bdf2*vnnn(1,1));
const double crhs24 = N[2]*(bdf0*vn(2,1) + bdf1*vnn(2,1) + bdf2*vnnn(2,1));
const double crhs25 = N[3]*(bdf0*vn(3,1) + bdf1*vnn(3,1) + bdf2*vnnn(3,1));
const double crhs26 = crhs22 + crhs23 + crhs24 + crhs25;
const double crhs27 = N[0]*(bdf0*vfrac(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*vfrac(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*vfrac(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*vfrac(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
const double crhs28 = crhs6*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1));
const double crhs29 = crhs8*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1));
const double crhs30 = crhs10*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1));
const double crhs31 = crhs28 + crhs29 + crhs30;
const double crhs32 = crhs13*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1)) + crhs14*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1)) + crhs15*(DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1));
const double crhs33 = -crhs22 - crhs23 - crhs24 - crhs25 + crhs27 - crhs28 - crhs29 - crhs30 + crhs32;
const double crhs34 = N[0]*(bdf0*vn(0,2) + bdf1*vnn(0,2) + bdf2*vnnn(0,2));
const double crhs35 = N[1]*(bdf0*vn(1,2) + bdf1*vnn(1,2) + bdf2*vnnn(1,2));
const double crhs36 = N[2]*(bdf0*vn(2,2) + bdf1*vnn(2,2) + bdf2*vnnn(2,2));
const double crhs37 = N[3]*(bdf0*vn(3,2) + bdf1*vnn(3,2) + bdf2*vnnn(3,2));
const double crhs38 = crhs34 + crhs35 + crhs36 + crhs37;
const double crhs39 = N[0]*(bdf0*vfrac(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*vfrac(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*vfrac(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*vfrac(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
const double crhs40 = crhs6*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2));
const double crhs41 = crhs8*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2));
const double crhs42 = crhs10*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2));
const double crhs43 = crhs40 + crhs41 + crhs42;
const double crhs44 = crhs13*(DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2)) + crhs14*(DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2)) + crhs15*(DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2));
const double crhs45 = -crhs34 - crhs35 - crhs36 - crhs37 + crhs39 - crhs40 - crhs41 - crhs42 + crhs44;
const double crhs46 = N[1]*crhs19;
const double crhs47 = crhs18*(DN(1,0)*crhs13 + DN(1,1)*crhs14 + DN(1,2)*crhs15);
const double crhs48 = N[2]*crhs19;
const double crhs49 = crhs18*(DN(2,0)*crhs13 + DN(2,1)*crhs14 + DN(2,2)*crhs15);
const double crhs50 = N[3]*crhs19;
const double crhs51 = crhs18*(DN(3,0)*crhs13 + DN(3,1)*crhs14 + DN(3,2)*crhs15);
rhs[0]=N[0]*crhs12 - N[0]*crhs16 + N[0]*crhs4 - N[0]*crhs5 - crhs17*crhs20 - crhs17*crhs21;
rhs[1]=N[0]*crhs26 - N[0]*crhs27 + N[0]*crhs31 - N[0]*crhs32 - crhs20*crhs33 - crhs21*crhs33;
rhs[2]=N[0]*crhs38 - N[0]*crhs39 + N[0]*crhs43 - N[0]*crhs44 - crhs20*crhs45 - crhs21*crhs45;
rhs[3]=N[1]*crhs12 - N[1]*crhs16 + N[1]*crhs4 - N[1]*crhs5 - crhs17*crhs46 - crhs17*crhs47;
rhs[4]=N[1]*crhs26 - N[1]*crhs27 + N[1]*crhs31 - N[1]*crhs32 - crhs33*crhs46 - crhs33*crhs47;
rhs[5]=N[1]*crhs38 - N[1]*crhs39 + N[1]*crhs43 - N[1]*crhs44 - crhs45*crhs46 - crhs45*crhs47;
rhs[6]=N[2]*crhs12 - N[2]*crhs16 + N[2]*crhs4 - N[2]*crhs5 - crhs17*crhs48 - crhs17*crhs49;
rhs[7]=N[2]*crhs26 - N[2]*crhs27 + N[2]*crhs31 - N[2]*crhs32 - crhs33*crhs48 - crhs33*crhs49;
rhs[8]=N[2]*crhs38 - N[2]*crhs39 + N[2]*crhs43 - N[2]*crhs44 - crhs45*crhs48 - crhs45*crhs49;
rhs[9]=N[3]*crhs12 - N[3]*crhs16 + N[3]*crhs4 - N[3]*crhs5 - crhs17*crhs50 - crhs17*crhs51;
rhs[10]=N[3]*crhs26 - N[3]*crhs27 + N[3]*crhs31 - N[3]*crhs32 - crhs33*crhs50 - crhs33*crhs51;
rhs[11]=N[3]*crhs38 - N[3]*crhs39 + N[3]*crhs43 - N[3]*crhs44 - crhs45*crhs50 - crhs45*crhs51;


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
