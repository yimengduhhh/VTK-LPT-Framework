#ifndef vtkLagrangianFinallyModel_h
#define vtkLagrangianFinallyModel_h

#include "vtkFiltersFlowPathsModule.h" // For export macro
#include "vtkLagrangianBasicIntegrationModel.h"

class VTKFILTERSFLOWPATHS_EXPORT vtkLagrangianFinallyModel
  : public vtkLagrangianBasicIntegrationModel
{
public:
  vtkTypeMacro(vtkLagrangianFinallyModel, vtkLagrangianBasicIntegrationModel);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkLagrangianFinallyModel* New();

  // Needed for multiple signatures polymorphism
  using Superclass::FunctionValues;

  /**
   * Evaluate the integration model velocity field
   * f at position x, using data from cell in dataSet with index cellId
   */
  int FunctionValues(vtkLagrangianParticle* particle, vtkDataSet* dataSet, vtkIdType cellId,
    double* weights, double* x, double* f) override;

protected:
  vtkLagrangianFinallyModel();
  ~vtkLagrangianFinallyModel() override;

  static double ShearReynolds(
    const double vorticity[3], double flowDensity, double particleDiameter, double dynVisc);

  static double ParticleReynolds(const double* flowVelocity, const double* particleVelocity,
    double dynVisc, double particleDiameter, double flowDensity);

  static double relativeVelocityX(const double* flowVelocity, const double* particleVelocity);
  static double relativeVelocityY(const double* flowVelocity, const double* particleVelocity);
  static double relativeVelocityZ(const double* flowVelocity, const double* particleVelocity);

  static double GetDragCoefficient(const double* flowVelocity, const double* particleVelocity,
    double dynVisc, double particleDiameter, double flowDensity);

  static double RelativeRotation(
    double Vorticity_x, double Vorticity_y, double Vorticity_z, double wpx, double wpy, double wpz);

  static double RotationReynolds(
    double flowDensity, double particleDiameter, double dynVisc, double RelativeRotation);

  static double GetAccelerationnumberAC(const double* flowVelocity, const double* particleVelocity,
    double particleDiameter, double modeltimestep, std::vector<double> LastVelCha);

  static double GetTorqueCoefficient(double Rerotation);

  static double GetCunninghamCoefficient(const double* flowVelocity, const double* particleVelocity,
    double dynVisc, double particleDiameter, double flowDensity);

  void GetDiffDiffVel(const double* flowVelocity, const double* particleVelocity,
    double* diffdiffvel, std::vector<double> LastVelCha);

  static double GetMagRelativeVelocity(const double* flowVelocity, const double* particleVelocity);

  static double GetDragNearWall(double distance, double particleDiameter, double Rep);

  static double GetShearNearWall(double distance, double particleDiameter, double Rep);

  static double Gaussian(double mean, double variance);

private:
  vtkLagrangianFinallyModel(const vtkLagrangianFinallyModel&) = delete;
  void operator=(const vtkLagrangianFinallyModel&) = delete;
};

#endif
