#include "vtkLagrangianFinallyModel.h"

#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkLagrangianParticle.h"
#include "vtkLagrangianParticleTracker.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStringArray.h"

#include <cstring>
#include <math.h>

#include <fstream>
#include <sstream>

vtkObjectFactoryNewMacro(vtkLagrangianFinallyModel);

vtkLagrangianFinallyModel::vtkLagrangianFinallyModel()
{

  this->SeedArrayNames->InsertNextValue("ParticleDiameter");
  this->SeedArrayComps->InsertNextValue(1);
  this->SeedArrayTypes->InsertNextValue(VTK_DOUBLE);
  this->SeedArrayNames->InsertNextValue("ParticleDensity");
  this->SeedArrayComps->InsertNextValue(1);
  this->SeedArrayTypes->InsertNextValue(VTK_DOUBLE);

  this->NumFuncs = 18;
  this->NumIndepVars = 19;
}

vtkLagrangianFinallyModel::~vtkLagrangianFinallyModel() = default;

void vtkLagrangianFinallyModel::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

int vtkLagrangianFinallyModel::FunctionValues(vtkLagrangianParticle* particle, vtkDataSet* dataSet,
  vtkIdType cellId, double* weights, double* x, double* f)
{

  std::fill(f, f + this->NumFuncs, 0.0);

  if (!particle)
  {
    vtkErrorMacro(<< "No particle to integrate");
    return 0;
  }

  if (!dataSet || cellId == -1)
  {
    vtkErrorMacro(<< "No cell or dataset to integrate the particle on. Dataset: " << dataSet
                  << " CellId:" << cellId);
    return 0;
  }

  double flowVelocity[3];
  if (this->GetFlowOrSurfaceDataNumberOfComponents(3, dataSet) != 3 ||
    !this->GetFlowOrSurfaceData(particle, 3, dataSet, cellId, weights, flowVelocity))
  {
    vtkErrorMacro(<< "Flow velocity is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Matida equations");
    return 0;
  }

  double flowDensity;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(4, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 4, dataSet, cellId, weights, &flowDensity))
  {
    vtkErrorMacro(<< "Flow density is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Matida equations");
    return 0;
  }

  double flowDynamicViscosity;
  if (this->IfManualViscosity == true)
  {
    flowDynamicViscosity = this->Viscosity;
  }
  else
  {

    if (this->GetFlowOrSurfaceDataNumberOfComponents(5, dataSet) != 1 ||
      !this->GetFlowOrSurfaceData(particle, 5, dataSet, cellId, weights, &flowDynamicViscosity))
    {
      vtkErrorMacro(<< "Flow dynamic viscosity is not set in source flow dataset or "
                       "has incorrect number of components, cannot use Matida equations");
      return 0;
    }
  }

  vtkDataArray* particleDiameters = vtkDataArray::SafeDownCast(this->GetSeedArray(6, particle));
  if (!particleDiameters)
  {
    vtkErrorMacro(<< "Particle diameter is not set in particle data, "
                     "cannot use Matida equations");
    return 0;
  }
  if (particleDiameters->GetNumberOfComponents() != 1)
  {
    vtkErrorMacro(<< "Particle diameter does not have the right number of components, "
                     "cannot use Matida equations");
    return 0;
  }
  double particleDiameter;
  particleDiameters->GetTuple(particle->GetSeedArrayTupleIndex(), &particleDiameter);

  vtkDataArray* particleDensities = vtkDataArray::SafeDownCast(this->GetSeedArray(7, particle));
  if (!particleDensities)
  {
    vtkErrorMacro(<< "Particle density is not set in particle data, "
                     "cannot use Matida equations");
    return 0;
  }
  if (particleDensities->GetNumberOfComponents() != 1)
  {
    vtkErrorMacro(<< "Particle densities does not have the right number of components, "
                     "cannot use Matida equations");
    return 0;
  }
  double particleDensity;
  particleDensities->GetTuple(particle->GetSeedArrayTupleIndex(), &particleDensity);

  double dxVelocitydy;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(8, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 8, dataSet, cellId, weights, &dxVelocitydy))
  {
    vtkErrorMacro(<< " dxVelocitydy is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dxVelocitydz;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(9, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 9, dataSet, cellId, weights, &dxVelocitydz))
  {
    vtkErrorMacro(<< "dxVelocitydz is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dyVelocitydx;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(10, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 10, dataSet, cellId, weights, &dyVelocitydx))
  {
    vtkErrorMacro(<< "dyVelocitydx is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dyVelocitydz;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(11, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 11, dataSet, cellId, weights, &dyVelocitydz))
  {
    vtkErrorMacro(<< "dyVelocitydz is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dzVelocitydx;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(12, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 12, dataSet, cellId, weights, &dzVelocitydx))
  {
    vtkErrorMacro(<< "dzVelocitydx is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dzVelocitydy;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(13, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 13, dataSet, cellId, weights, &dzVelocitydy))
  {
    vtkErrorMacro(<< "dzVelocitydy is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double Temperature;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(14, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 14, dataSet, cellId, weights, &Temperature))
  {
    vtkErrorMacro(<< "Temperature is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }
  if (this->IfManualTemperature == true)
  {
    Temperature = this->ModelTemperature;
  }

  double dpressuredx;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(15, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 15, dataSet, cellId, weights, &dpressuredx))
  {
    vtkErrorMacro(<< "dpressuredx is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dpressuredy;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(16, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 16, dataSet, cellId, weights, &dpressuredy))
  {
    vtkErrorMacro(<< "dpressuredy is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double dpressuredz;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(17, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 17, dataSet, cellId, weights, &dpressuredz))
  {
    vtkErrorMacro(<< "dpressuredz is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double turbulentKE;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(18, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 18, dataSet, cellId, weights, &turbulentKE))
  {
    vtkErrorMacro(<< "turbulentKE is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double turbulentDR;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(19, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 19, dataSet, cellId, weights, &turbulentDR))
  {
    vtkErrorMacro(<< "turbulentKE is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double fillId;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(20, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 20, dataSet, cellId, weights, &fillId))
  {
    vtkErrorMacro(<< "fillId is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double RMSE_X;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(21, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 21, dataSet, cellId, weights, &RMSE_X))
  {
    vtkErrorMacro(<< "fillId is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double RMSE_Y;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(22, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 22, dataSet, cellId, weights, &RMSE_Y))
  {
    vtkErrorMacro(<< "fillId is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double RMSE_Z;
  if (this->GetFlowOrSurfaceDataNumberOfComponents(23, dataSet) != 1 ||
    !this->GetFlowOrSurfaceData(particle, 23, dataSet, cellId, weights, &RMSE_Z))
  {
    vtkErrorMacro(<< "fillId is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Shear equations");
    return 0;
  }

  double Vorticity[3];
  if (this->GetFlowOrSurfaceDataNumberOfComponents(25, dataSet) != 3 ||
    !this->GetFlowOrSurfaceData(particle, 25, dataSet, cellId, weights, Vorticity))
  {
    vtkErrorMacro(<< "Vorticity is not set in source flow dataset or "
                     "has incorrect number of components, cannot use Matida equations");
    return 0;
  }

  Vorticity[0] = dzVelocitydy - dyVelocitydz;
  Vorticity[1] = dxVelocitydz - dzVelocitydx;
  Vorticity[2] = dyVelocitydx - dxVelocitydy;

  f[3] = 0;
  f[4] = 0;
  f[5] = 0;

  double mass;
  mass = 3.1415926 * particleDiameter * particleDiameter * particleDiameter * particleDensity / 6.0;

  if (this->IfOpenTurbulentDispersion == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {
    turbulentKE = 0.5 * (RMSE_X * RMSE_X + RMSE_Y * RMSE_Y + RMSE_Z * RMSE_Z);
    turbulentDR = std::pow(2.0 / 3.0 * turbulentKE, 1.5) / 0.02;
    double C_T = 0.24;
    double C_L = 3.0;
    double sigma = sqrt(2.0 / 3.0 * turbulentKE);
    double T_L = C_T * pow(sigma, 2) / turbulentDR;
    double L_E = C_T * T_L * sigma;
    double dt = this->ModelTimeStep;
    double R_L = exp(-dt / T_L);
    double deltaR[3];
    deltaR[0] = (flowVelocity[0] - particle->GetVelocity()[0]) * dt;
    deltaR[1] = (flowVelocity[1] - particle->GetVelocity()[1]) * dt;
    deltaR[2] = (flowVelocity[2] - particle->GetVelocity()[2]) * dt;
    double magdeltaR = sqrt(deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2]);
    double fDeltaR = exp(-magdeltaR / L_E);
    double gDeltaR = fDeltaR * (1.0 - magdeltaR / 2.0 / L_E);
    double R_E[3];
    double R_P[3];

    for (int i = 0; i < 3; i++)
    {
      R_E[i] = (fDeltaR - gDeltaR) * deltaR[i] * deltaR[i] / magdeltaR / magdeltaR + gDeltaR;
      R_P[i] = R_L * R_E[i];
    }
    double Wiener[3];

    for (int i = 0; i < 3; i++)
    {

      Wiener[i] = vtkLagrangianFinallyModel::Gaussian(0, dt);
    }
    for (int i = 0; i < 3; i++)
    {
      f[15 + i] = R_P[i] * x[15 + i] + sigma * sqrt(1 - R_P[i] * R_P[i]) * Wiener[i];

      flowVelocity[i] = flowVelocity[i] + f[15 + i];
    }
  }

  double relativeVelocityxx =
    vtkLagrangianFinallyModel::relativeVelocityX(flowVelocity, particle->GetVelocity());
  double relativeVelocityyy =
    vtkLagrangianFinallyModel::relativeVelocityY(flowVelocity, particle->GetVelocity());
  double relativeVelocityzz =
    vtkLagrangianFinallyModel::relativeVelocityZ(flowVelocity, particle->GetVelocity());

  double forcegravityx = 0;
  double forcegravityy = 0;
  double forcegravityz = 0;
  if (this->IfOpenGravity == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    forcegravityx = 0;
    forcegravityy = 0;
    forcegravityz = 9.8;

    f[3] = f[3] + forcegravityx;
    f[4] = f[4] + forcegravityy;
    f[5] = f[5] + forcegravityz;
  }
  else if (this->IfOpenGravityx == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    forcegravityx = -9.8;
    forcegravityy = 0;
    forcegravityz = 0;

    f[3] = f[3] + forcegravityx;
    f[4] = f[4] + forcegravityy;
    f[5] = f[5] + forcegravityz;
  }
  else if (this->IfOpenGravityz == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    forcegravityx = 0;
    forcegravityy = 0;
    forcegravityz = -9.8;

    f[3] = f[3] + forcegravityx;
    f[4] = f[4] + forcegravityy;
    f[5] = f[5] + forcegravityz;
  }
  else
  {
    forcegravityx = 0;
    forcegravityy = 0;
    forcegravityz = 0;
  }

  double forcebuoyancyx = 0;
  double forcebuoyancyy = 0;
  double forcebuoyancyz = 0;
  if (this->IfOpenBuoyancy == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    if (this->IfOpenGravityx == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
    {
      forcebuoyancyx = 0;
      forcebuoyancyy = 0;
      forcebuoyancyz = -flowDensity / particleDensity * 9.8;
    }
    else if (this->IfOpenGravityz == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
    {
      forcebuoyancyx = 0;
      forcebuoyancyy = 0;
      forcebuoyancyz = flowDensity / particleDensity * 9.8;
    }
    else
    {
      forcebuoyancyx = 0;
      forcebuoyancyy = flowDensity / particleDensity * 9.8;
      forcebuoyancyz = 0;
    }

    f[3] = f[3] + forcebuoyancyx;
    f[4] = f[4] + forcebuoyancyy;
    f[5] = f[5] + forcebuoyancyz;
  }

  double forcedragx = 0;
  double forcedragy = 0;
  double forcedragz = 0;
  if (this->IfOpenDrag == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    double drag = vtkLagrangianFinallyModel::GetDragCoefficient(
      flowVelocity, particle->GetVelocity(), flowDynamicViscosity, particleDiameter, flowDensity);
    double MagRelativeVel =
      vtkLagrangianFinallyModel::GetMagRelativeVelocity(flowVelocity, particle->GetVelocity());
    forcedragx = 3.0 * flowDensity * drag * relativeVelocityxx * MagRelativeVel / 4.0 /
      particleDensity / particleDiameter;
    forcedragy = 3.0 * flowDensity * drag * relativeVelocityyy * MagRelativeVel / 4.0 /
      particleDensity / particleDiameter;
    forcedragz = 3.0 * flowDensity * drag * relativeVelocityzz * MagRelativeVel / 4.0 /
      particleDensity / particleDiameter;

    f[3] = f[3] + forcedragx;
    f[4] = f[4] + forcedragy;
    f[5] = f[5] + forcedragz;
  }

  double forcecunningx = 0;
  double forcecunningy = 0;
  double forcecunningz = 0;
  if (this->IfOpenCunninghamdrag == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {
    double drag = vtkLagrangianFinallyModel::GetDragCoefficient(
      flowVelocity, particle->GetVelocity(), flowDynamicViscosity, particleDiameter, flowDensity);
    double cunningcoefficient = vtkLagrangianFinallyModel::GetCunninghamCoefficient(
      flowVelocity, particle->GetVelocity(), flowDynamicViscosity, particleDiameter, flowDensity);
    double MagRelativeVel =
      vtkLagrangianFinallyModel::GetMagRelativeVelocity(flowVelocity, particle->GetVelocity());
    forcecunningx = 3.0 * flowDensity * drag / cunningcoefficient * relativeVelocityxx *
      MagRelativeVel / 4.0 / particleDensity / particleDiameter;
    forcecunningy = 3.0 * flowDensity * drag / cunningcoefficient * relativeVelocityyy *
      MagRelativeVel / 4.0 / particleDensity / particleDiameter;
    forcecunningz = 3.0 * flowDensity * drag / cunningcoefficient * relativeVelocityzz *
      MagRelativeVel / 4.0 / particleDensity / particleDiameter;
    f[3] = f[3] + forcecunningx;
    f[4] = f[4] + forcecunningy;
    f[5] = f[5] + forcecunningz;
  }

  double forceshearliftx = 0;
  double forceshearlifty = 0;
  double forceshearliftz = 0;
  if (this->IfOpenShearlift == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    double PReynolds;
    double SReynolds;
    double CorrectionFunction;
    PReynolds = vtkLagrangianFinallyModel::ParticleReynolds(
      flowVelocity, particle->GetVelocity(), flowDynamicViscosity, particleDiameter, flowDensity);
    SReynolds = vtkLagrangianFinallyModel::ShearReynolds(
      Vorticity, flowDensity, particleDiameter, flowDynamicViscosity);
    double beta = 0.5 * SReynolds / PReynolds;
    if (PReynolds < 40)
    {
      CorrectionFunction =
        (1 - 0.3314 * std::sqrt(beta)) * std::exp(-PReynolds / 10.0) + 0.3314 * std::sqrt(beta);
    }
    else
    {
      CorrectionFunction = 0.0524 * std::sqrt(beta * PReynolds);
    }

    double LiftCoefficient = 4.1126 * CorrectionFunction / std::sqrt(SReynolds);

    double ForceCoefficient;
    ForceCoefficient = 3.0 * flowDensity * LiftCoefficient / 4.0 / particleDensity;

    double RelativeVel[3];
    for (int i = 0; i < 3; i++)
    {
      RelativeVel[i] = flowVelocity[i] - particle->GetVelocity()[i];
    }

    double VelCrossVor[3];
    vtkMath::Cross(RelativeVel, Vorticity, VelCrossVor);
    forceshearliftx = ForceCoefficient * VelCrossVor[0];
    forceshearlifty = ForceCoefficient * VelCrossVor[1];
    forceshearliftz = ForceCoefficient * VelCrossVor[2];

    f[3] = f[3] + forceshearliftx;
    f[4] = f[4] + forceshearlifty;
    f[5] = f[5] + forceshearliftz;
  }

  double relative_omega[3];
  relative_omega[0] = Vorticity[0] / 2.0 - x[6];
  relative_omega[1] = Vorticity[1] / 2.0 - x[7];
  relative_omega[2] = Vorticity[2] / 2.0 - x[8];

  double mag_omega;
  mag_omega = vtkLagrangianFinallyModel::RelativeRotation(
    Vorticity[0], Vorticity[1], Vorticity[2], x[6], x[7], x[8]);

  double Rerotation;
  Rerotation = vtkLagrangianFinallyModel::RotationReynolds(
    flowDensity, particleDiameter, flowDynamicViscosity, mag_omega);

  double Clr = vtkLagrangianFinallyModel::GetTorqueCoefficient(Rerotation);

  f[6] =
    15.0 / (16 * 3.1415926) * (flowDensity / particleDensity) * Clr * mag_omega * relative_omega[0];
  f[7] =
    15.0 / (16 * 3.1415926) * (flowDensity / particleDensity) * Clr * mag_omega * relative_omega[1];
  f[8] =
    15.0 / (16 * 3.1415926) * (flowDensity / particleDensity) * Clr * mag_omega * relative_omega[2];

  double forcerotationx = 0;
  double forcerotationy = 0;
  double forcerotationz = 0;
  if (this->IfOpenRotationlift == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    double Reparticle;
    Reparticle = vtkLagrangianFinallyModel::ParticleReynolds(
      flowVelocity, particle->GetVelocity(), flowDynamicViscosity, particleDiameter, flowDensity);

    double LiftRotationCoefficient;
    LiftRotationCoefficient = 0.45 +
      (Rerotation / Reparticle - 0.45) *
        std::exp(-0.05684 * pow(Rerotation, 0.4) * pow(Reparticle, 0.3));

    double CrelativeVelocity[3];
    for (int i = 0; i < 3; i++)
    {
      CrelativeVelocity[i] = particle->GetVelocity()[i] - flowVelocity[i];
    }
    double CrelativeSpeed = vtkMath::Norm(CrelativeVelocity);
    double TotalCoefficient;
    TotalCoefficient = 3 * flowDensity * LiftRotationCoefficient * CrelativeSpeed / 4.0 /
      particleDiameter / mag_omega / particleDensity;

    forcerotationx = TotalCoefficient *
      (relative_omega[1] * relativeVelocityzz - relative_omega[2] * relativeVelocityyy);
    forcerotationy = TotalCoefficient *
      (relative_omega[2] * relativeVelocityxx - relative_omega[0] * relativeVelocityzz);
    forcerotationz = TotalCoefficient *
      (relative_omega[0] * relativeVelocityyy - relative_omega[1] * relativeVelocityxx);

    f[3] = f[3] + forcerotationx;
    f[4] = f[4] + forcerotationy;
    f[5] = f[5] + forcerotationz;
  }

  double forcepressuregradx = 0;
  double forcepressuregrady = 0;
  double forcepressuregradz = 0;
  if (this->IfOpenPressuregrad == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {

    forcepressuregradx = -dpressuredx / particleDensity;
    forcepressuregrady = -dpressuredy / particleDensity;
    forcepressuregradz = -dpressuredz / particleDensity;

    f[3] = f[3] + forcepressuregradx;
    f[4] = f[4] + forcepressuregrady;
    f[5] = f[5] + forcepressuregradz;
  }

  double addmassaccx = 0;
  double addmassaccy = 0;
  double addmassaccz = 0;
  if (this->IfOpenAddmass == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {
    if (particle->GetIfFirstStep() == true)
    {
      particle->SetIfFirstStep(false);
    }
    else
    {
      std::vector<double> LastVelCha = particle->GetLastVelocityCha();
      double accelerationnumberAC;
      accelerationnumberAC = vtkLagrangianFinallyModel::GetAccelerationnumberAC(
        flowVelocity, particle->GetVelocity(), particleDiameter, this->ModelTimeStep, LastVelCha);
      double coefficientA;
      coefficientA = 2.1 - 0.132 / (accelerationnumberAC * accelerationnumberAC + 0.12);
      double* diffdiffvel = new double[3];
      this->GetDiffDiffVel(flowVelocity, particle->GetVelocity(), diffdiffvel, LastVelCha);
      addmassaccx =
        0.5 * coefficientA * flowDensity / particleDensity * diffdiffvel[0] / this->ModelTimeStep;
      addmassaccy =
        0.5 * coefficientA * flowDensity / particleDensity * diffdiffvel[1] / this->ModelTimeStep;
      addmassaccz =
        0.5 * coefficientA * flowDensity / particleDensity * diffdiffvel[2] / this->ModelTimeStep;

      f[3] = f[3] + addmassaccx;
      f[4] = f[4] + addmassaccy;
      f[5] = f[5] + addmassaccz;
      delete[] diffdiffvel;
    }
  }

  double brownianx = 0;
  double browniany = 0;
  double brownianz = 0;
  if (this->IfOpenBrownian == vtkLagrangianBasicIntegrationModel::FORCE_OPEN)
  {
    double randomnumx = vtkMath::Random();
    double randomnumy = vtkMath::Random();
    double randomnumz = vtkMath::Random();
    double Ccun = vtkLagrangianFinallyModel::GetCunninghamCoefficient(
      flowVelocity, particle->GetVelocity(), flowDynamicViscosity, particleDiameter, flowDensity);
    double kBoltz = 1.3806488e-23;
    double brownianCoeff1 = sqrt(216 * kBoltz / 3.1415926 / (particleDensity / flowDensity) /
      (particleDensity / flowDensity) / Ccun);
    double brownianCoeff2 = sqrt(flowDynamicViscosity / flowDensity / flowDensity * Temperature /
      pow(particleDiameter, 5) / this->ModelTimeStep);
    brownianx = brownianx + randomnumx * brownianCoeff1 * brownianCoeff2;
    browniany = browniany + randomnumy * brownianCoeff1 * brownianCoeff2;
    brownianz = brownianz + randomnumz * brownianCoeff1 * brownianCoeff2;

    f[3] = f[3] + brownianx;
    f[4] = f[4] + browniany;
    f[5] = f[5] + brownianz;
  }

  for (int i = 0; i < 3; i++)
  {
    f[i] = x[i + 3];
  }

  f[9] = flowVelocity[0];
  f[10] = flowVelocity[1];
  f[11] = flowVelocity[2];

  f[12] = flowVelocity[0] - particle->GetVelocity()[0];
  f[13] = flowVelocity[1] - particle->GetVelocity()[1];
  f[14] = flowVelocity[2] - particle->GetVelocity()[2];

  std::vector<double> VelCha(3);
  VelCha.clear();
  VelCha.push_back(flowVelocity[0] - particle->GetVelocity()[0]);
  VelCha.push_back(flowVelocity[1] - particle->GetVelocity()[1]);
  VelCha.push_back(flowVelocity[2] - particle->GetVelocity()[2]);
  particle->SetLastVelocityCha(VelCha);

  return 1;
}

double vtkLagrangianFinallyModel::ShearReynolds(
  const double vorticity[3], double flowDensity, double particleDiameter, double dynVisc)
{
  double magVor;
  magVor = vtkMath::Norm(vorticity);
  return flowDensity * particleDiameter * particleDiameter * magVor / dynVisc;
}

double vtkLagrangianFinallyModel::ParticleReynolds(const double* flowVelocity,
  const double* particleVelocity, double dynVisc, double particleDiameter, double flowDensity)
{
  if (dynVisc == 0)
  {
    return -1.0 * std::numeric_limits<double>::infinity();
  }
  double relativeVelocity[3];
  for (int i = 0; i < 3; i++)
  {
    relativeVelocity[i] = flowVelocity[i] - particleVelocity[i];
  }
  double relativeSpeed = vtkMath::Norm(relativeVelocity);
  double reynolds = flowDensity * relativeSpeed * particleDiameter / dynVisc;
  return reynolds;
}

double vtkLagrangianFinallyModel::relativeVelocityX(
  const double* flowVelocity, const double* particleVelocity)
{

  return (flowVelocity[0] - particleVelocity[0]);
}

double vtkLagrangianFinallyModel::relativeVelocityY(
  const double* flowVelocity, const double* particleVelocity)
{

  return (flowVelocity[1] - particleVelocity[1]);
}

double vtkLagrangianFinallyModel::relativeVelocityZ(
  const double* flowVelocity, const double* particleVelocity)
{

  return (flowVelocity[2] - particleVelocity[2]);
}

double vtkLagrangianFinallyModel::GetDragCoefficient(const double* flowVelocity,
  const double* particleVelocity, double dynVisc, double particleDiameter, double flowDensity)
{
  if (dynVisc == 0)
  {
    return -1.0 * std::numeric_limits<double>::infinity();
  }
  double relativeVelocity[3];
  for (int i = 0; i < 3; i++)
  {
    relativeVelocity[i] = flowVelocity[i] - particleVelocity[i];
  }
  double relativeSpeed = vtkMath::Norm(relativeVelocity);
  double reynolds = flowDensity * relativeSpeed * particleDiameter / dynVisc;

  if (reynolds == 0)
  {
    return 0;
  }

  if (reynolds < 0.5)
  {
    return 24.0 / reynolds;
  }
  else if (reynolds >= 1000)
  {
    return 0.44;
  }
  else
  {
    return 24.0 * (1 + 0.15 * std::pow(reynolds, 0.687)) / reynolds;
  }
}

double vtkLagrangianFinallyModel::RelativeRotation(
  double Vorticity_x, double Vorticity_y, double Vorticity_z, double wpx, double wpy, double wpz)
{
  double xRotation = Vorticity_x / 2 - wpx;
  double yRotation = Vorticity_y / 2 - wpy;
  double zRotation = Vorticity_z / 2 - wpz;
  return std::sqrt(xRotation * xRotation + yRotation * yRotation + zRotation * zRotation);
}

double vtkLagrangianFinallyModel::RotationReynolds(
  double flowDensity, double particleDiameter, double dynVisc, double RelativeRotation)
{
  return flowDensity * particleDiameter * particleDiameter * RelativeRotation / dynVisc;
}

double vtkLagrangianFinallyModel::GetAccelerationnumberAC(const double* flowVelocity,
  const double* particleVelocity, double particleDiameter, double modeltimestep,
  std::vector<double> LastVelCha)
{
  double a1 = (flowVelocity[0] - particleVelocity[0]) * (flowVelocity[0] - particleVelocity[0]);
  double a2 = (flowVelocity[1] - particleVelocity[1]) * (flowVelocity[1] - particleVelocity[1]);
  double a3 = (flowVelocity[2] - particleVelocity[2]) * (flowVelocity[2] - particleVelocity[2]);
  double newdiffvel[3];
  for (int i = 0; i < 3; i++)
  {
    newdiffvel[i] = flowVelocity[i] - particleVelocity[i];
  }
  double olddiffvel[3];

  for (int i = 0; i < 3; i++)
  {
    olddiffvel[i] = LastVelCha[i];
  }

  double cha[3];
  for (int i = 0; i < 3; i++)
  {
    cha[i] = newdiffvel[i] - olddiffvel[i];
  }

  double chaMag = sqrt(cha[0] * cha[0] + cha[1] * cha[1] + cha[2] * cha[2]);
  return (a1 + a2 + a3) / particleDiameter / chaMag * modeltimestep;
}

double vtkLagrangianFinallyModel::GetTorqueCoefficient(double Rerotation)
{
  double TorqueCoefficient;
  if (Rerotation < 32)
  {
    return 64 * 3.1415926 / Rerotation;
  }
  else
  {
    TorqueCoefficient = 12.9 / std::sqrt(Rerotation) + 128.4 / Rerotation;
    return TorqueCoefficient;
  }
}

double vtkLagrangianFinallyModel::GetCunninghamCoefficient(const double* flowVelocity,
  const double* particleVelocity, double dynVisc, double particleDiameter, double flowDensity)
{
  double meanfreepath = 0.00000007;
  double a1 = 0.4 * std::exp(-1.1 * particleDiameter / 2.0 / meanfreepath);
  return (1 + 2 * meanfreepath / particleDiameter * (1.257 + a1));
}

void vtkLagrangianFinallyModel::GetDiffDiffVel(const double* flowVelocity,
  const double* particleVelocity, double* diffdiffvel, std::vector<double> LastVelCha)
{
  double newdiffvel[3];
  for (int i = 0; i < 3; i++)
  {
    newdiffvel[i] = flowVelocity[i] - particleVelocity[i];
  }
  double olddiffvel[3];
  for (int i = 0; i < 3; i++)
  {
    olddiffvel[i] = LastVelCha[i];
  }

  for (int i = 0; i < 3; i++)
  {
    diffdiffvel[i] = newdiffvel[i] - olddiffvel[i];
  }
}

double vtkLagrangianFinallyModel::GetMagRelativeVelocity(
  const double* flowVelocity, const double* particleVelocity)
{
  double relativeVel[3];
  for (int i = 0; i < 3; i++)
  {
    relativeVel[i] = flowVelocity[i] - particleVelocity[i];
  }
  double b = vtkMath::Norm(relativeVel);
  return b;
}

double vtkLagrangianFinallyModel::Gaussian(double mean, double variance)
{

  double v1, v2, s;
  double x;

  do
  {
    double u1 = vtkMath::Random();
    double u2 = vtkMath::Random();
    v1 = 2 * u1 - 1;
    v2 = 2 * u2 - 1;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1 || s == 0);
  x = v1 * sqrt(-2 * log(s) / s);
  return x * variance + mean;
}
