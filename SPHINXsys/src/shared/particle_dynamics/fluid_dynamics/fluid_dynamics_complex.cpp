/**
 * @file 	fluid_dynamics_complex.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		void ShearThinningViscousAccelerationWithWall::Interaction(size_t index_i, Real dt)
		{
			ShearThinningViscousAccelerationInner::Interaction(index_i, dt);

			Real rho_i = this->rho_n_[index_i];
			const Vecd &vel_i = this->vel_n_[index_i];
			Real mu_i = this->mu_shear_[index_i]; 

			Vecd acceleration(0), vel_derivative(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(this->wall_Vol_[k]);
				StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
				Neighborhood &contact_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					acceleration += 2.0 * mu_i * vel_derivative * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] / rho_i;
					//acceleration += 2.0 * this->mu_ * vel_derivative * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] / rho_i;
				}
			}

			this->dvel_dt_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(BaseBodyRelationInner &inner_relation,
											   BaseBodyRelationContact &contact_relation)
			: ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(inner_relation, contact_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation)
			: TransportVelocityCorrectionComplex(complex_relation.inner_relation_,
												 complex_relation.contact_relation_) {}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation,
											   BaseBodyRelationContact &extra_contact_relation)
			: ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(complex_relation, extra_contact_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::prepareContactData()
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::Interaction(size_t index_i, Real dt)
		{
			TransportVelocityCorrectionInner::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];

			Vecd acceleration_trans(0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * Vol_k[index_j] * nablaW_ij / rho_i;
				}
			}

			/** correcting particle position */
			if (surface_indicator_[index_i] == 0)
				pos_n_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxation<PressureRelaxationInnerOldroyd_B>::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];
					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * nablaW_ij * Vol_k[index_j] / rho_i;
				}
			}

			dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		void DensityRelaxationWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxation<DensityRelaxationInnerOldroyd_B>::Interaction(index_i, dt);

			Vecd vel_i = vel_n_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
				StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];

					Matd velocity_gradient = -SimTK::outer((vel_i - vel_ave_k[index_j]), nablaW_ij) * Vol_k[index_j] * 2.0;
					stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
								   (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
				}
			}
			dtau_dt_[index_i] += stress_rate;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//