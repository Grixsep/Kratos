//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_BOUNDARY_ROTATION_UTILITY
#define KRATOS_MPM_BOUNDARY_ROTATION_UTILITY

// system includes

// external includes

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "utilities/coordinate_transformation_utilities.h"

// Application includes
#include "particle_mechanics_application_variables.h"

namespace Kratos {

///@addtogroup ParticleMechanicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/* A utility to rotate the local contributions of certain nodes to the system matrix,
which is required to apply slip conditions (roller-type support) in arbitrary directions to the boundary nodes.*/
template<class TLocalMatrixType, class TLocalVectorType>
class MPMBoundaryRotationUtility: public CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double> {
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of MPMBoundaryRotationUtility
	KRATOS_CLASS_POINTER_DEFINITION(MPMBoundaryRotationUtility);

	using CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>::Rotate;

	typedef Node NodeType;

	typedef Geometry< Node > GeometryType;

	///@}
	///@name Life Cycle
	///@{

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Displacement DOFs are assumed to be the first mDomainSize rows in each block of rows.
	 * @param rVariable Kratos variable used to flag nodes where local system contributions will be rotated. All nodes with rVariable != Zero will be rotated.
	 */
	MPMBoundaryRotationUtility(
        const unsigned int DomainSize,
		const unsigned int BlockSize,
		const Variable<double>& rVariable):
    CoordinateTransformationUtils<TLocalMatrixType,TLocalVectorType,double>(DomainSize,BlockSize,SLIP), mrFlagVariable(rVariable)
	{}

	/// Destructor.
	~MPMBoundaryRotationUtility() override {}

	/// Assignment operator.
	MPMBoundaryRotationUtility& operator=(MPMBoundaryRotationUtility const& rOther) {}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	/// Rotate the local system contributions so that they are oriented with each node's normal.
	/**
	 @param rLocalMatrix Local system matrix
	 @param rLocalVector Local RHS vector
	 @param rGeometry A reference to the element's (or condition's) geometry
	 */
	void Rotate(
        TLocalMatrixType& rLocalMatrix,
		TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const override
	{
		if (this->GetBlockSize() == this->GetDomainSize()) // irreducible case
		{
			if (this->GetDomainSize() == 2) this->template RotateAuxPure<2>(rLocalMatrix,rLocalVector,rGeometry);
			else if (this->GetDomainSize() == 3) this->template RotateAuxPure<3>(rLocalMatrix,rLocalVector,rGeometry);
		}
		else // mixed formulation case
		{
			if (this->GetDomainSize() == 2) this->template RotateAux<2,3>(rLocalMatrix,rLocalVector,rGeometry);
			else if (this->GetDomainSize() == 3) this->template RotateAux<3,4>(rLocalMatrix,rLocalVector,rGeometry);
		}

	}

	/// RHS only version of Rotate
	void RotateRHS(
        TLocalVectorType& rLocalVector,
		GeometryType& rGeometry) const
	{
		this->Rotate(rLocalVector,rGeometry);
	}

    // Auxiliary function to clear friction-related flags --
    // MUST be called before (re-)building the RHS in a given non-linear iteration
    static void ClearFrictionFlag(const ModelPart &rModelPart) {
        KRATOS_TRY
        // Loop over the grid nodes performed to clear INLET flag for indicating that nodal friction has alr been set
        // and remove the accumulated normal forces
        for(NodeType &curr_node : rModelPart.Nodes()){
            curr_node.SetLock();
            curr_node.Reset(INLET);
            curr_node.FastGetSolutionStepValue(MPM_NORMAL_FORCE) = 0.0;
            curr_node.UnSetLock();
        }
        KRATOS_CATCH( "" )
    }

	/// Apply roler type boundary conditions to the rotated local contributions.
	/** This function takes the rotated local system contributions so each
	 node's displacement are expressed using a base oriented with its normal
	 and imposes that the normal displacement is equal to the mesh displacement in
	 the normal direction.
	 */
	void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const override
	{
		const unsigned int LocalSize = rLocalVector.size();

		if (LocalSize > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				if(this->IsSlip(rGeometry[itNode]) )
				{
					// We fix the first displacement dof (normal component) for each rotated block
					unsigned int j = itNode * this->GetBlockSize();

					// Get the displacement of the boundary mesh, this does not assume that the mesh is moving.
					// If the mesh is moving, need to consider the displacement of the moving mesh into account.
					const array_1d<double,3> & displacement = rGeometry[itNode].FastGetSolutionStepValue(DISPLACEMENT);

					// Get Normal Vector of the boundary
					array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);

                    // Zero out row/column corresponding to normal displacement DoF except diagonal term (set to 1)
                    // Applied IFF the local matrix passed is not empty [otherwise does nothing -- RHS only case]
                    if (rLocalMatrix.size1() != 0) {
                        for( unsigned int i = 0; i < LocalSize; ++i)
                        {
                            rLocalMatrix(i,j) = 0.0;
                            rLocalMatrix(j,i) = 0.0;
                        }
                        rLocalMatrix(j, j) = 1.0; // set diagonal term to 1.0
                    }

                    /// Computation of nodal reaction forces due to conforming SLIP BC -- use MPM_NORMAL_FORCE
                    /// as frame of reference is aligned with node normals [& not global frame of ref per REACTION]
                    // Accumulate RHS values along normal direction to FORCE_RESIDUAL [nodal reaction forces due to conforming SLIP]
                    // -- when converged, RHS ~= 0 -> FORCE_RESIDUAL = -RHS value [computed without adding reaction force]

                    rGeometry[itNode].SetLock();
                    rGeometry[itNode].FastGetSolutionStepValue(MPM_NORMAL_FORCE) -= rLocalVector[j];
                    rGeometry[itNode].UnSetLock();

                    // Set value of normal displacement at node directly to the normal displacement of the boundary mesh
					rLocalVector[j] = inner_prod(rN,displacement);

                    // TODO: refactor? [e.g. set INLET flag and local flag, unlock node and do another if based on local flag]
                    // Prescribe a constant force of FRICTION_FORCE Newtons in the 1st tangential direction to all nodes
                    // INLET flag (reset in MPMResidualBasedBossakScheme.InitializeSolutionStep) ensures that friction
                    // is applied exactly once per node
                    rGeometry[itNode].SetLock();
                    if(!rGeometry[itNode].Is(INLET)) {
                        // obtain nodal velocity and rotate it to the same frame of reference as the local geometry
                        array_1d<double, 3> nodal_velocity = ZeroVector(3);

                        nodal_velocity = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);

                        // TODO: make threshold configurable?
                        this->RotateAndNormalizeVector(nodal_velocity, rGeometry[itNode], VELOCITY_THRESHOLD);

                        // apply friction force in opposite direction of tangential velocity components
                        for (unsigned dim = 1; dim < this->GetDomainSize(); dim++) {
                            rLocalVector[j + dim] -= FRICTION_FORCE * nodal_velocity[dim];
                        }

                        rGeometry[itNode].Set(INLET);
                    }
                    rGeometry[itNode].UnSetLock();
				}
			}
		}
	}

	/// RHS only version of ApplySlipCondition
	void ApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const override
	{
        // creates an empty dummy matrix to pass into the 'full' ApplySlipCondition -- this dummy matrix is
        // ignored, effectively only updating the RHS
        TLocalMatrixType dummyMatrix;
        this->ApplySlipCondition(dummyMatrix, rLocalVector, rGeometry);
	}

	// An extra function to distinguish the application of slip in element considering penalty imposition
	void ElementApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		// If it is not a penalty element, do as standard
		// Otherwise, if it is a penalty element, don't do anything
		if (!this->IsPenalty(rGeometry))
		{
			this->ApplySlipCondition(rLocalMatrix, rLocalVector, rGeometry);
		}
	}

	// An extra function to distinguish the application of slip in element considering penalty imposition (RHS Version)
	void ElementApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		// If it is not a penalty element, do as standard
		// Otherwise, if it is a penalty element, don't do anything
		if (!this->IsPenalty(rGeometry))
		{
			this->ApplySlipCondition(rLocalVector, rGeometry);
		}
	}

	// An extra function to distinguish the application of slip in condition considering penalty imposition
	void ConditionApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		// If it is not a penalty condition, do as standard
		if (!this->IsPenalty(rGeometry))
		{
			this->ApplySlipCondition(rLocalMatrix, rLocalVector, rGeometry);
		}
		// Otherwise, do the following modification
		else
		{
			const unsigned int LocalSize = rLocalVector.size();

			if (LocalSize > 0)
			{
				const unsigned int block_size = this->GetBlockSize();
				TLocalMatrixType temp_matrix = ZeroMatrix(rLocalMatrix.size1(),rLocalMatrix.size2());
				for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
				{
					if(this->IsSlip(rGeometry[itNode]) )
					{
						// We fix the first displacement dof (normal component) for each rotated block
						unsigned int j = itNode * block_size;

						// Copy all normal value in LHS to the temp_matrix
                        // [ does nothing for dummy rLocalMatrix (size1() == 0) -- RHS only case ]
						for (unsigned int i = j; i < rLocalMatrix.size1(); i+= block_size)
						{
							temp_matrix(i,j) = rLocalMatrix(i,j);
							temp_matrix(j,i) = rLocalMatrix(j,i);
						}

						// Remove all other value in RHS than the normal component
						for(unsigned int i = j; i < (j + block_size); ++i)
						{
							if (i!=j) rLocalVector[i] = 0.0;
						}

                        // Add 'friction' force //

                        // Only apply to nodes currently containing material points
                        if (rGeometry[itNode].FastGetSolutionStepValue(NODAL_MASS) >= std::numeric_limits<double>::epsilon()) {
                            // obtain nodal velocity and rotate it to the same frame of reference as the local geometry
                            // [ since mpc contributions are accumulated at nodes, friction force is directed against
                            //    nodal velocity value                                                                ]

//                            KRATOS_WATCH(rGeometry[itNode].GetInitialPosition());

                            array_1d<double, 3> nodal_velocity = ZeroVector(3);

                            nodal_velocity = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);

                            this->RotateAndNormalizeVector(nodal_velocity, rGeometry[itNode], VELOCITY_THRESHOLD);

                            // TODO: do this in a condition object (makes more sense since later you'll need the nodal force)
                            // obtain friction contribution of at boundary particle [currently fixed FRICTION_FORCE] and extrapolate to nodes
                            Vector shape_fn = row(rGeometry.ShapeFunctionsValues(), 0);
                            double nodal_friction_contribution = shape_fn[itNode] * FRICTION_FORCE;

//                                                    KRATOS_WATCH(j);
//                                                    KRATOS_WATCH(nodal_velocity);

                            // apply friction force in opposite direction of tangential velocity components
                            for (unsigned dim = 1; dim < this->GetDomainSize(); dim++) {
                                rLocalVector[j + dim] -= nodal_friction_contribution * nodal_velocity[dim];

                                //                            KRATOS_WATCH(rLocalVector[j+dim]);
                            }
                        }
                    }
				}
                // All entries in penalty matrix zeroed out except for normal component
                // [ no effect in case of empty dummy rLocalMatrix ]
				rLocalMatrix = temp_matrix;
			}
		}
	}

	// An extra function to distinguish the application of slip in condition considering penalty imposition (RHS Version)
	void ConditionApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
        // creates an empty dummy matrix to pass into the 'full' ConditionApplySlipCondition -- this dummy matrix is
        // ignored, effectively only updating the RHS
        TLocalMatrixType dummyMatrix;
        this->ConditionApplySlipCondition(dummyMatrix, rLocalVector, rGeometry);
	}

	// Checking whether it is normal element or penalty element
	bool IsPenalty(GeometryType& rGeometry) const
	{
		bool is_penalty = false;
		for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
		{
			if(this->IsSlip(rGeometry[itNode]) )
			{
				const double identifier = rGeometry[itNode].FastGetSolutionStepValue(mrFlagVariable);
				const double tolerance  = 1.e-6;
				if (identifier > 1.00 + tolerance)
				{
					is_penalty = true;
					break;
				}
			}
		}

		return is_penalty;
	}

	/// Same functionalities as RotateVelocities, just to have a clear function naming
	virtual	void RotateDisplacements(ModelPart& rModelPart) const
	{
		this->RotateVelocities(rModelPart);
	}

	/// Transform nodal displacement to the rotated coordinates (aligned with each node's normal)
	/// The name is kept to be Rotate Velocities, since it is currently a derived class of coordinate_transformation_utilities in the core
	void RotateVelocities(ModelPart& rModelPart) const override
	{
		TLocalVectorType displacement(this->GetDomainSize());
		TLocalVectorType Tmp(this->GetDomainSize());

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
		#pragma omp parallel for firstprivate(displacement,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				//this->RotationOperator<TLocalMatrixType>(Rotation,);
				if(this->GetDomainSize() == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 3; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(rRot,displacement);
					for(unsigned int i = 0; i < 3; i++) rDisplacement[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 2; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(rRot,displacement);
					for(unsigned int i = 0; i < 2; i++) rDisplacement[i] = Tmp[i];
				}
			}
		}
	}

	/// Same functionalities as RecoverVelocities, just to have a clear function naming
	virtual void RecoverDisplacements(ModelPart& rModelPart) const
	{
		this->RecoverVelocities(rModelPart);
	}

	/// Transform nodal displacement from the rotated system to the original configuration
	/// The name is kept to be Recover Velocities, since it is currently a derived class of coordinate_transformation_utilities in the core
	void RecoverVelocities(ModelPart& rModelPart) const override
	{
		TLocalVectorType displacement(this->GetDomainSize());
		TLocalVectorType Tmp(this->GetDomainSize());

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
		#pragma omp parallel for firstprivate(displacement,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				if(this->GetDomainSize() == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 3; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(trans(rRot),displacement);
					for(unsigned int i = 0; i < 3; i++) rDisplacement[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					this->LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
					for(unsigned int i = 0; i < 2; i++) displacement[i] = rDisplacement[i];
					noalias(Tmp) = prod(trans(rRot),displacement);
					for(unsigned int i = 0; i < 2; i++) rDisplacement[i] = Tmp[i];
				}
			}
		}
	}

    // Constraints velocities and accelerations normal to SLIP condition to zero
    void ConstraintDerivatives(array_1d<double, 3>& rCurrentVelocity,
                               array_1d<double, 3>& rCurrentAcceleration,
                               const NodeType& rNode){
        // rotate to local normal-tangential frame
        this->RotateVector(rCurrentVelocity, rNode);
        this->RotateVector(rCurrentAcceleration, rNode);

        rCurrentVelocity[0] = 0;
        rCurrentAcceleration[0] = 0;

        // rotate back to global coordinates
        this->RotateVector(rCurrentVelocity, rNode, true);
        this->RotateVector(rCurrentAcceleration, rNode, true);
    }

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const override
	{
		std::stringstream buffer;
		buffer << "MPMBoundaryRotationUtility";
		return buffer.str();
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override
	{
		rOStream << "MPMBoundaryRotationUtility";
	}

	/// Print object's data.
	void PrintData(std::ostream& rOStream) const override {}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	const Variable<double>& mrFlagVariable;

    const double FRICTION_FORCE = 0000;
    const double VELOCITY_THRESHOLD = 1e-10;

	///@}
	///@name Member Variables
	///@{

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{


    /// Helper function to rotate a 3-vector to and from the coordinate system defined by the NORMAL defined at rNode
    /**
     @param rVector Vector to be rotated
     @param rNode A reference to the node associated with the vector
     @param toGlobalCoordinates If true, instead rotates the vector back to the global coordinates [default: false]
     */
    void RotateVector( array_1d<double, 3> &rVector,
                       const Node &rNode,
                       const bool toGlobalCoordinates = false,
                       const double threshold = std::numeric_limits<double>::epsilon() ) const {
        array_1d<double, 3> rotated_nodal_vector = ZeroVector(3);
        BoundedMatrix<double, 3, 3> rotation_matrix = ZeroMatrix(3);

        this->LocalRotationOperatorPure(rotation_matrix, rNode);
        noalias(rotated_nodal_vector) = prod(toGlobalCoordinates ? trans(rotation_matrix) : rotation_matrix, rVector);

        rVector = rotated_nodal_vector;
    }

    /// Additionally normalizes the rotated vector
    /**
     @param rVector Vector to be rotated
     @param rNode A reference to the node associated with the vector
     @param toGlobalCoordinates If true, instead rotates the vector back to the global coordinates [default: false]
     @param threshold Value below which the value of the component is considered 0 [ default value: machine epsilon ]
     */
    void RotateAndNormalizeVector( array_1d<double, 3> &rVector,
                      const Node &rNode,
                      const bool toGlobalCoordinates = false,
                      const double threshold = std::numeric_limits<double>::epsilon() ) const {
        RotateVector(rVector, rNode, toGlobalCoordinates);

        // Check if velocity is close to zero [ALL components below threshold]
        bool is_zero_vector = true;

        // Evaluates to false if ANY of the rVector component is above threshold
        for (auto component : rVector) {
            is_zero_vector = is_zero_vector && (abs(component) < threshold);
        }

        // If not close to zero, rotate the vector and obtain its norm
        // Otherwise do nothing
        if(!is_zero_vector){
            this->Normalize(rVector);
        }
    }
	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::istream& operator >>(std::istream& rIStream,
		MPMBoundaryRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType>
inline std::ostream& operator <<(std::ostream& rOStream,
		const MPMBoundaryRotationUtility<TLocalMatrixType, TLocalVectorType>& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_MPM_BOUNDARY_ROTATION_UTILITY
