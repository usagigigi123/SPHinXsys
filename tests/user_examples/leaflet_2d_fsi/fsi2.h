/**
* @file 	fsi2.h
* @brief 	This is the case file for the test of fluid - structure interaction.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#ifndef FSI2_CASE_H
#define FSI2_CASE_H

#include "sphinxsys.h"
#include "linear.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.001;
Real DH = 20. * scale;                  /**< Channel height. */
Real R = DH;                            /**< Sinus radius. */
Real DL1 = 2. * DH;					/**< Channel upstream length. */
Real DL2 = 4.5 * DH;					/**< Channel downstream length. */	
Real W = 120. * scale;                  /**< Leaflet width. */
Real resolution_ref = DH / 60.0;				/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;			/**< Boundary width, determined by specific layer of boundary particles. */
Vec2d insert_circle_center(DL1 + R, DH); 
/** Beam related parameters. */
Vec2d beam_fixed_point(DL1, DH);
Real cos_a = sqrt(2.)/2.;
Real sin_a = sqrt(2.)/2.;  
Real bl = 26.0 * scale;	                /**< Length of the beam. */           
Real bh = 1.0 * scale;	                /**< Thickness of the beam. */
Real refinement_ratio = resolution_ref / (bh / 6.);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL1 + DL2 + BW, DH + R + BW));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 890.0;											  /**< Density. */
Real U_f = 23.44e-3/60./DH/W;												  /**< Characteristic velocity. */
Real U_max = 2. * U_f;
Real c_f = 10.0 * U_max;										  /**< Speed of sound. */
Real mu_f = 4.3e-3; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1000.0; /**< Reference density.*/
Real poisson = 0.49; /**< Poisson ratio.*/
Real Youngs_modulus = 1.5e6 * std::pow(0.3e-3/bh, 3);
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, DH));
	water_block_shape.push_back(Vecd(DL1 + DL2, DH));
	water_block_shape.push_back(Vecd(DL1 + DL2, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

	return water_block_shape;
}
/** create a beam shape */
Real hbh = bh / 2.0;
Vec2d BLB(beam_fixed_point[0] - hbh * sin_a, beam_fixed_point[1] - hbh * cos_a);
Vec2d BLT(beam_fixed_point[0] + hbh * sin_a, beam_fixed_point[1] + hbh * cos_a);
Vec2d BRB(beam_fixed_point[0] + bl * cos_a - hbh * sin_a, beam_fixed_point[1] - bl * sin_a - hbh * cos_a);
Vec2d BRT(beam_fixed_point[0] + bl * cos_a + hbh * sin_a, beam_fixed_point[1] - bl * sin_a + hbh * cos_a);
//Beam observer location
StdVec<Vecd> beam_observation_location = {0.5 * (BRT + BRB)};
std::vector<Vecd> createBeamShape()
{
	std::vector<Vecd> beam_shape;
	beam_shape.push_back(BLB);
	beam_shape.push_back(BLT);
	beam_shape.push_back(BRT);
	beam_shape.push_back(BRB);
	beam_shape.push_back(BLB);

	return beam_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL1 + DL2 + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL1 + DL2 + BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL1 + DL2 + 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL1 + DL2 + 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

	return inner_wall_shape;
}
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
//----------------------------------------------------------------------
//	Define case dependent geometrices
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon_.addACircle(insert_circle_center, R, 100, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
	}
};
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addACircle(insert_circle_center, R + BW, 100, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
		multi_polygon_.addACircle(insert_circle_center, R, 100, ShapeBooleanOps::sub);
	}
};
class Insert : public MultiPolygonShape
{
public:
	explicit Insert(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createBeamShape(), ShapeBooleanOps::add);
	}
};
/** create the beam base as constrain shape. */
std::vector<Vecd> createConstraintShape()
{
	std::vector<Vecd> constraint_shape;
	constraint_shape.push_back(BLB);
	constraint_shape.push_back(BLT);
	constraint_shape.push_back(Vecd(beam_fixed_point[0] + BW * cos_a + hbh * sin_a, beam_fixed_point[1] - BW * sin_a + hbh * cos_a));
	constraint_shape.push_back(Vecd(beam_fixed_point[0] + BW * cos_a - hbh * sin_a, beam_fixed_point[1] - BW * sin_a - hbh * cos_a));
	constraint_shape.push_back(BLB);

	return constraint_shape;
}
MultiPolygon createBeamBaseShape()
{
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(createConstraintShape(), ShapeBooleanOps::add);
	return multi_polygon;
}

/** Case dependent inflow boundary condition. */
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody &fluid_body, BodyAlignedBoxByCell &aligned_box_part)
		: InflowBoundaryCondition(fluid_body, aligned_box_part),
		  u_ave_(0), u_ref_(U_f), t_ref(1.0) {}
	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0)
		{
			u = 1.5 * u_ave_ * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
			//u = 1.5 * u_ref_ * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
			v = 0.0;
		}
		return Vecd(u, v);
	}
	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		run_time = run_time - std::floor(run_time / t_ref) * t_ref;
		Real Q = interpolate(run_time);
		u_ave_ = Q * 1.e-3/60./DH/W;
	}
};
/** fluid observer particle generator */
class FluidObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit FluidObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** A line of measuring points at the entrance of the channel. */
		Real x1 = DL1 - 1.5 * DH;
		Real x2 = DL1 + 4 * DH;
		size_t number_observation_points = 21;
		/** the measuring locations */
		for (size_t i = 0; i < number_observation_points; ++i)
		{
			positions_.push_back(Vec2d(x1, DH * (Real)i / (Real)(number_observation_points - 1)));
			positions_.push_back(Vec2d(x2, DH * (Real)i / (Real)(number_observation_points - 1)));
		}
	}
};

#endif //FSI2_CASE_H
