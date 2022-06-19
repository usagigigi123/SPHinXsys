#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real scale = 1.0;
Real PL = 25.0*scale;									  /** Length of the square plate. */
Real PH = 25.0*scale;									  /** Width of the square plate. */
Real PT = 1.0*scale;									  /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 1.0, 0.0);				  /** Pseudo-normal. */						  /** Particle number in the direction of the length */
Real resolution_ref = 1.0*scale; /** Initial reference particle spacing. */
//int BWD = 4;
//Real BW = resolution_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
Real level_set_refinement_ratio = resolution_ref / (0.05 * PT);
/** Domain bounds of the system. */
TriangleMeshShapeSTL triangle_mesh_geometry_shape("./input/leaflet_1mm_1.stl", Vecd(0., 0., 0.)*scale, scale);
//TriangleMeshShapeSTL triangle_mesh_geometry_shape("./input/aorta_valve_center_1_0.5.stl", Vecd(0., 0., -3.6358), scale);
//TriangleMeshShapeSTL triangle_mesh_geometry_shape("./input/aorta_valve_center_2.stl", Vecd(0., 0., -PH/2.0)*scale, scale);
BoundingBox system_domain_bounds = triangle_mesh_geometry_shape.findBounds();

/** For material properties of the solid. */
Real rho0_s = 1100.0e-6;				   /** Normalized density. */
Real Youngs_modulus = 1.62e4; /** Normalized Youngs Modulus. */
Real poisson = 0.0;				   /** Poisson ratio. */
Real physical_viscosity = 200.0;   /** physical damping, here we choose the same value as numerical viscosity. */

Real q = 0; /** Total distributed load. */
Real time_to_full_external_force = 1.0;

Real gravitational_acceleration = 1.0 * 1.0e-2;

/** Define application dependent particle generator for thin structure. */
class ShellBodyFromMesh : public ThinStructure
{
public:
    ShellBodyFromMesh(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape, SharedPtr<SPHAdaptation> particle_adaptation)
	: ThinStructure(system, body_name, particle_adaptation)
    {
        body_shape_.add<LevelSetShape>(this, triangle_mesh_shape, true, false);
        // set the body domain bounds because it is not set by default
        BoundingBox bounds = body_shape_.findBounds();
        setBodyDomainBounds(bounds);
    }
};

/**
 * application dependent initial condition
 */
class PlateDynamicsInitialCondition
	: public thin_structure_dynamics::ShellDynamicsInitialCondition
{
public:
	explicit PlateDynamicsInitialCondition(SolidBody &plate)
		: thin_structure_dynamics::ShellDynamicsInitialCondition(plate){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial pseudo-normal. */
		n_0_[index_i] = n_0;
		//n_0_[index_i] = Vec3d(pos_0_[index_i][0] / radius_mid_surface, 0.0, pos_0_[index_i][2] / radius_mid_surface);
		n_[index_i] = n_0;
		pseudo_n_[index_i] = n_0_[index_i];
		transformation_matrix_[index_i] = getTransformationMatrix(n_0_[index_i]);
	};
};

/** Define the boundary geometry. */
/*
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (type == "cylinder")
		{   if (base_particles_->pos_n_[index_i][1] < BW)
            //if (base_particles_->pos_n_[index_i][2] < radius_mid_surface * sin(-17.5 / 180.0 * Pi))
		    {
			    body_part_particles_.push_back(index_i);
		    }
		}
		if (type == "circle")
		{
			Real radius = PH/2.0;
		    Real a = sqrt(base_particles_->pos_n_[index_i][0]*base_particles_->pos_n_[index_i][0]
                     +base_particles_->pos_n_[index_i][2]*base_particles_->pos_n_[index_i][2]);
            if (a >= radius - BW)
            {
                body_part_particles_.push_back(index_i);
            }
		}
		if (type == "plate")
		{
			if (base_particles_->pos_n_[index_i][0] < -PL/2. + BW || base_particles_->pos_n_[index_i][2] < -PH/2. + BW ||
			    base_particles_->pos_n_[index_i][0] > PL/2.-BW || base_particles_->pos_n_[index_i][2] > PH/2.-BW)
		    {
			    body_part_particles_.push_back(index_i);
		    }
		}
	};
};
/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
public:
	explicit TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		//return global_acceleration_;
		//return global_acceleration_* 0.5 * (1.0 + sin(Pi * current_time / time_to_full_external_force - 0.5 * Pi));
		return current_time < time_to_full_external_force
				   ? current_time * global_acceleration_ / time_to_full_external_force
				   : global_acceleration_;
	}
};

/** Define application dependent observer particle generator. */
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		//positions_volumes_.push_back(std::make_pair(Vecd(0.5 * PL, 0.5 * PH, 0.0), 0.0));
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	//SPHSystem system(system_domain_bounds, particle_spacing_ref);
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** create a Cylinder body. */
	//ThinStructure plate_body(system, "PlateBody", makeShared<SPHAdaptation>(1.15, 1.0));
	/** create particles for the elastic body. */
	//ShellParticles plate_body_particles(plate_body,
	//									   makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
	//									   makeShared<PlateParticleGenerator>(), thickness);
	
	
	shared_ptr<SPHAdaptation> shell_adaptation = makeShared<SPHAdaptation>(1.15, 1.0, 0.75, level_set_refinement_ratio);
	ShellBodyFromMesh plate_body (system, "Shell_plate", triangle_mesh_geometry_shape, shell_adaptation);
	ShellParticles plate_body_particles(
		plate_body,
		makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
		//makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson),
		makeShared<ShellParticleGeneratorLattice>(PT), PT
	);

	/** Define Observer. */
	ObserverBody plate_observer(system, "PlateObserver");
	ObserverParticles observer_particles(plate_observer, makeShared<ObserverParticleGenerator>());

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner plate_body_inner(plate_body);
	BodyRelationContact plate_observer_contact(plate_observer, {&plate_body});

	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);


	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(plate_body);
	// A  Physics relaxation step. 
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(plate_body_inner, PT, level_set_refinement_ratio);
	
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	write_states.writeToFile(0.0);
	plate_body.updateCellLinkedList();
	
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
			//write_states.writeToFile(ite_p);
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	relaxation_step_inner.mid_surface_bounding_.calculateNormalDirection(); 

	// update pos_0_
	std::copy(
		plate_body.base_particles_->pos_n_.begin(),
		plate_body.base_particles_->pos_n_.end(),
		plate_body.base_particles_->pos_0_.begin()
	);

	write_states.writeToFile(ite_p);
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
    
	/** Common particle dynamics. */
	//TimeDependentExternalForce external_force(Vec3d(0.0, 0.0, q / (PT * rho0_s) - gravitational_acceleration));
	
	TimeDependentExternalForce external_force(Vec3d(0.0, gravitational_acceleration, 0.0));
	//TimeDependentExternalForce external_force(Vec3d(0.0,  0.0, gravitational_acceleration));
	TimeStepInitialization initialize_external_force(plate_body, external_force);

	
	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** initial condition */
	PlateDynamicsInitialCondition plate_initial_pseudo_normal(plate_body);
	/** Corrected configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration(plate_body_inner);
	/** Time step size calculation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(plate_body);
	/** active-passive stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(plate_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(plate_body_inner);
	/** Constrain the Boundary. */
	//BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	//solid_dynamics::ConstrainSolidBodyRegion constrain_holder(plate_body, boundary_geometry);
	
	TriangleMeshShapeSTL holder_shape("./input/leaflet_1mm_base.stl", Vecd(0.0, 0.0, 0.0)*scale, scale);
	BodyRegionByParticle boundary_geometry(plate_body, "Holder", holder_shape);
	thin_structure_dynamics::ConstrainShellBodyRegion constrain_holder(plate_body, boundary_geometry);
	
	
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		plate_position_damping(plate_body_inner, 0.5, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		plate_rotation_damping(plate_body_inner, 0.5, "AngularVelocity", physical_viscosity);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	plate_initial_pseudo_normal.parallel_exec();
	corrected_configuration.parallel_exec();

	/**
	* From here the time stepping begins.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);

	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 2.0;
	//Real end_time = 2.0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			plate_position_damping.parallel_exec(dt);
			plate_rotation_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

			plate_body.updateCellLinkedList();
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;


	return 0;
}
