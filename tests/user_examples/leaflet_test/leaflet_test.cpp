/**
 * @file 	test_3d_shell_particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate shell particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation by levelset for a complex thin structures geometry (3D).
 * @author 	Dong Wu and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_geometry = "./input/leaflet_1mm_1.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 1.;								 
Real PH = 25.0*scale;									  												
Vec3d domain_lower_bound(-PH, 0.0, -PH);
Vec3d domain_upper_bound(PH, 22.0*scale, PH);
Real resolution_ref = 1.*scale; 
Real PT = 1.*scale;
// level set resolution much higher than that of particles is required
Real level_set_refinement_ratio = resolution_ref / (0.1 * PT);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ImportedShellModel: public ComplexShape
{
public:
	explicit ImportedShellModel(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_geometry, Vecd(0)*scale, 1.0*scale);
	}
};

//For material properties of the solid.
Real rho0_s = 1100.0;				   
Real Youngs_modulus = 1.6e4;
Real poisson = 0.0;				   
Real physical_viscosity = 200.0;   

Real q = 0; 
Real time_to_full_external_force = 0.5;

Real gravitational_acceleration = 100.0 * 1.0e0;

class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name), body_shape_(body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
		
		body_shape_.add<TriangleMeshShapeSTL>("./input/leaflet_1mm_base.stl", Vecd(0.0, 0.0, 0.0)*scale, scale);
	};
	virtual ~BoundaryGeometry(){};

private:
    ComplexShape body_shape_;
	void tagManually(size_t index_i)
	{
		if (body_shape_.checkContain(base_particles_->pos_n_[index_i]))
		{
			body_part_particles_.push_back(index_i);
		}
	};
};


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


int main()
{
	//----------------------------------------------------------------------
	//	Build up a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	InOutput in_output(system);
    //----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody plate_body(system, makeShared<ImportedShellModel>("ImportedShellModel"));
	plate_body.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(plate_body);
	//here dummy linear elastic solid is use because no solid dynamics in particle relaxation
	plate_body.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	plate_body.generateParticles<ThickSurfaceParticleGeneratorLattice>(PT);
	//plate_body.addBodyStateForRecording<Vecd>("NormalDirection");

    //----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	//MeshRecordingToPlt write_mesh_cell_linked_list(in_output, imported_model, imported_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner plate_body_inner(plate_body);

	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizeParticlePosition  random_imported_model_particles(plate_body);
	// A  Physics relaxation step. 
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(plate_body_inner, PT, level_set_refinement_ratio);
	relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(plate_body_inner, PT);
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
	//relaxation_step_inner.mid_surface_bounding_.calculateNormalDirection(); 
    shell_normal_prediction.exec();
	// update pos_0_
//	std::copy(
//		plate_body.base_particles_->pos_n_.begin(),
//		plate_body.base_particles_->pos_n_.end(),
//		plate_body.base_particles_->pos_0_.begin()
//	);

	write_states.writeToFile(ite_p);
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
    
	//TimeDependentExternalForce external_force(Vec3d(0.0, 0.0, q / (PT * rho0_s) - gravitational_acceleration));
	TimeDependentExternalForce external_force(Vec3d(0.0, gravitational_acceleration, 0.0));
	//TimeDependentExternalForce external_force(Vec3d(0.0,  0.0, gravitational_acceleration));
	TimeStepInitialization initialize_external_force(plate_body, external_force);

	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration(plate_body_inner);
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(plate_body);
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(plate_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(plate_body_inner);
	//BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	//solid_dynamics::ConstrainSolidBodyRegion constrain_holder(plate_body, boundary_geometry);
	
	//TriangleMeshShapeSTL holder_shape("./input/leaflet_1mm_base.stl", Vecd(0.0, 0.0, 0.0)*scale, scale);
	//BodyRegionByParticle boundary_geometry(plate_body, "Holder", holder_shape);
	BodyRegionByParticle boundary_geometry(plate_body, makeShared<TriangleMeshShapeSTL>("./input/leaflet_1mm_base.stl", Vecd(0.0, 0.0, 0.0)*scale, scale));

	//BoundaryGeometry boundary_geometry(plate_body, "Holder");
	thin_structure_dynamics::ConstrainShellBodyRegion constrain_holder(plate_body, boundary_geometry);
	
	
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		plate_position_damping(0.5, plate_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		plate_rotation_damping(0.5, plate_body_inner, "AngularVelocity", physical_viscosity);

	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	//plate_initial_pseudo_normal.parallel_exec();
	corrected_configuration.parallel_exec();

	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);

	int ite = 0;
	Real end_time = 1.0;
	//Real end_time = 2.0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

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
			dt = 0.01*computing_time_step_size.parallel_exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

			//plate_body.updateCellLinkedList();
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