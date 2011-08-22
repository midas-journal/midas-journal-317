/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkDMFFDLabeledPointSetRegistrationFilter.txx,v $
Language:  C++

Date:      $Date: $
Version:   $Revision: $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDMFFDLabeledPointSetRegistrationFilter_txx
#define __itkDMFFDLabeledPointSetRegistrationFilter_txx

#include "itkDMFFDLabeledPointSetRegistrationFilter.h"

#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkJensenHavrdaCharvatTsallisLabeledPointSetMetric.h"

#include "vnl/vnl_math.h"

namespace itk {

/**
 * N3BiasFieldScaleCostFunction class definitions
 */
template<class TFilter>
DMFFDLabeledPointSetRegistrationCostFunction<TFilter>
::DMFFDLabeledPointSetRegistrationCostFunction()
{
  this->m_Filter = NULL;
}

template<class TFilter>
DMFFDLabeledPointSetRegistrationCostFunction<TFilter>
::~DMFFDLabeledPointSetRegistrationCostFunction()
{
}

template<class TFilter>
typename DMFFDLabeledPointSetRegistrationCostFunction<TFilter>::MeasureType
DMFFDLabeledPointSetRegistrationCostFunction<TFilter>
::GetValue( const ParametersType & parameters ) const
{
  // This function should not be needed as the optimizer
  // only calls GetValueAndDerivative()
  return 0;
}

template<class TFilter>
void
DMFFDLabeledPointSetRegistrationCostFunction<TFilter>
::GetDerivative( const ParametersType & parameters,
  DerivativeType & derivative ) const
{
  // This function should not be needed as the optimizer
  // only calls GetValueAndDerivative()
}

template<class TFilter>
void
DMFFDLabeledPointSetRegistrationCostFunction<TFilter>
::GetValueAndDerivative( const ParametersType & parameters, MeasureType & value,
  DerivativeType & derivative ) const
{
  /**
   * Warp the moving points
   */
  typename TFilter::MovingPointSetType::Pointer warpedPoints =
    TFilter::MovingPointSetType::New();
  warpedPoints->Initialize();

  typename TFilter::ControlPointLatticeType::Pointer controlPointLattice =
    TFilter::ControlPointLatticeType::New();
  controlPointLattice->SetOrigin( this->m_Filter->
    m_TotalDeformationFieldControlPointLattice->GetOrigin() );
  controlPointLattice->SetSpacing( this->m_Filter->
    m_TotalDeformationFieldControlPointLattice->GetSpacing() );
  controlPointLattice->SetRegions( this->m_Filter->
    m_TotalDeformationFieldControlPointLattice->GetLargestPossibleRegion() );
  controlPointLattice->SetDirection( this->m_Filter->
    m_TotalDeformationFieldControlPointLattice->GetDirection() );
  controlPointLattice->Allocate();

  unsigned long index = 0;
  ImageRegionIterator<typename TFilter::ControlPointLatticeType> It(
    controlPointLattice, controlPointLattice->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename TFilter::ControlPointLatticeType::PixelType vector;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      vector[d] = parameters[index++];
      }
    It.Set( vector );
    }

  typename TFilter::ControlPointFilterType::Pointer bspliner
    = TFilter::ControlPointFilterType::New();
  bspliner->SetSplineOrder( this->m_Filter->m_SplineOrder );
  bspliner->SetInput( controlPointLattice  );
  bspliner->SetOrigin( this->m_Filter->m_Origin );
  bspliner->SetSize( this->m_Filter->m_Size );
  bspliner->SetSpacing( this->m_Filter->m_Spacing );
  bspliner->SetDirection( this->m_Filter->m_Direction );

  unsigned long count = 0;
  for( unsigned int n = 0; n < this->m_Filter->GetInput( 1 )->
    GetNumberOfPoints(); n++ )
    {
    typename TFilter::MovingPointType inputPoint;
    this->m_Filter->GetInput( 1 )->GetPoint( n, &inputPoint );

    typename TFilter::ControlPointFilterType::PointType point;
    point.CastFrom( inputPoint );

    typename TFilter::VectorType vector;
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      continue;
      }
    inputPoint += vector;

    typename TFilter::MovingPointSetType::PixelType inputLabel = 1;
    this->m_Filter->GetInput( 1 )->GetPointData( n, &inputLabel );

    warpedPoints->SetPoint( count, inputPoint );
    warpedPoints->SetPointData( count, inputLabel );
    count++;
    }

  /**
   * Set up the point-set function
   */
  if( this->m_Filter->m_Verbose )
    {
    std::cout << "     *** Current annealing temperature: "
      << this->m_Filter->m_PointSetSigma *
         vcl_sqrt( vcl_pow( this->m_Filter->m_AnnealingRate, static_cast<RealType>(
         this->m_Filter->m_Optimizer->GetCurrentIteration() ) ) )
      << " ***" << std::endl;
    }


  typedef JensenHavrdaCharvatTsallisLabeledPointSetMetric
    <typename TFilter::FixedPointSetType> PointSetFunctionType;
  typename PointSetFunctionType::Pointer pointSetFunction
    = PointSetFunctionType::New();

  pointSetFunction->SetUseWithRespectToTheMovingPointSet( true );
  pointSetFunction->SetUseAnisotropicCovariances(
    this->m_Filter->m_UseAnisotropicCovariances );
  pointSetFunction->SetUseInputAsSamples(
    this->m_Filter->m_UseInputAsSamples );
  pointSetFunction->SetUseRegularizationTerm( this->m_Filter->
    m_UseRegularizationTerm );
  pointSetFunction->SetAlpha( this->m_Filter->m_Alpha );

  pointSetFunction->SetFixedPointSet( this->m_Filter->GetInput( 0 ) );
  pointSetFunction->SetFixedPointSetSigma( this->m_Filter->m_PointSetSigma *
    vcl_sqrt( vcl_pow( this->m_Filter->m_AnnealingRate, static_cast<RealType>(
    this->m_Filter->m_Optimizer->GetCurrentIteration() ) ) ) );
  pointSetFunction->SetFixedKernelSigma( this->m_Filter->m_KernelSigma );
  pointSetFunction->SetFixedCovarianceKNeighborhood(
    this->m_Filter->m_CovarianceKNeighborhood );
  pointSetFunction->SetFixedEvaluationKNeighborhood(
    this->m_Filter->m_EvaluationKNeighborhood );
  pointSetFunction->SetNumberOfFixedSamples( this->m_Filter->
    m_NumberOfFixedSamples );

  pointSetFunction->SetMovingPointSet( warpedPoints );
  pointSetFunction->SetMovingPointSetSigma( this->m_Filter->m_PointSetSigma *
    vcl_sqrt( vcl_pow( this->m_Filter->m_AnnealingRate, static_cast<RealType>(
    this->m_Filter->m_Optimizer->GetCurrentIteration() ) ) ) );
  pointSetFunction->SetMovingKernelSigma( this->m_Filter->m_KernelSigma );
  pointSetFunction->SetMovingCovarianceKNeighborhood(
    this->m_Filter->m_CovarianceKNeighborhood );
  pointSetFunction->SetMovingEvaluationKNeighborhood(
    this->m_Filter->m_EvaluationKNeighborhood );
  pointSetFunction->SetNumberOfMovingSamples( this->m_Filter->
    m_NumberOfMovingSamples );

  pointSetFunction->Initialize();

  typename PointSetFunctionType::DefaultTransformType::ParametersType
    zeroParameters;
  zeroParameters.Fill( 0 );
  typename PointSetFunctionType::DerivativeType pointSetGradient;
  typename PointSetFunctionType::MeasureType zeroValue;
  pointSetFunction->GetValueAndDerivative( zeroParameters, zeroValue,
    pointSetGradient );

  value = 0.0;
  for( unsigned int n = 0; n < zeroValue.Size(); n++ )
    {
    value += zeroValue[n];
    }

  typename TFilter::BSplinePointSetType::Pointer fieldPoints
    = TFilter::BSplinePointSetType::New();
  fieldPoints->Initialize();
  typename TFilter::BSplineFilterType::WeightsContainerType::Pointer weights
    = TFilter::WeightsContainerType::New();
  weights->Initialize();

  count = 0;
  for( unsigned int n = 0; n < pointSetFunction->GetNumberOfValues(); n++ )
    {
    typename TFilter::BSplinePointSetType::PointType fieldPoint;
    typename TFilter::VectorType gradient;

    typename TFilter::MovingPointType warpedPoint;
    warpedPoints->GetPoint( n, &warpedPoint );

    typename TFilter::WarpedPixelType label = 1;
    warpedPoints->GetPointData( n, &label );

    bool isInside = true;
    for( unsigned d = 0; d < TFilter::Dimension; d++ )
      {
      gradient[d] = pointSetGradient(n, d);
      fieldPoint[d] = warpedPoint[d];
      if( fieldPoint[d] <= this->m_Filter->m_Origin[d] ||
          fieldPoint[d] >= this->m_Filter->m_Origin[d] +
          this->m_Filter->m_Spacing[d] *
          static_cast<RealType>( this->m_Filter->m_Size[d] - 1 ) )
        {
        isInside = false;
        }
      }
    if( isInside )
      {
      fieldPoints->SetPoint( count, fieldPoint );
      fieldPoints->SetPointData( count, gradient );
      if( this->m_Filter->m_LabelWeights.size() == 0 )
        {
        weights->InsertElement( count, 1.0 );
        }
      else
        {
        weights->InsertElement( count, this->m_Filter->m_LabelWeights[label] );
        }
      count++;
      }
    }

//   itkDebugMacro( "Normalizing gradient values such that each gradient has" <<
//     "a norm equivalent to the voxel spacing." );
//
//   RealType sumSquaredNorm = 0.0;
//   for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
//     {
//     typename TFilter::VectorType gradient;
//     fieldPoints->GetPointData( i, &gradient );
//     sumSquaredNorm += gradient.GetSquaredNorm();
//     }
//   typename TFilter::VectorType V;
//   RealType sigma;
//   for( unsigned int i = 0; i < Dimension; i++ )
//     {
//     V[i] = this->m_Filter->m_Spacing[i];
//     }
//   if( Dimension == 2 )
//     {
//     sigma = V.GetNorm();
//     }
//   else if( Dimension == 3 )
//     {
//     sigma = V.GetNorm() / vcl_sqrt( 2.0 );
//     }
//
//   RealType gradientNormalizationFactor = sigma * vcl_sqrt
//     ( static_cast<RealType>( Dimension *
//     fieldPoints->GetNumberOfPoints() ) / sumSquaredNorm );
//
//   for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
//     {
//     typename TFilter::VectorType gradient;
//     fieldPoints->GetPointData( i, &gradient );
//     fieldPoints->SetPointData( i, gradient * gradientNormalizationFactor );
//     }

  typename TFilter::BSplineFilterType::ArrayType nlevels;
  typename TFilter::BSplineFilterType::ArrayType ncps;

  nlevels.Fill( 1 );
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    ncps[d] = this->m_Filter->m_TotalDeformationFieldControlPointLattice
      ->GetLargestPossibleRegion().GetSize()[d];
    }

  typename TFilter::BSplineFilterType::Pointer bsplineFilter =
    TFilter::BSplineFilterType::New();

  bsplineFilter->SetInput( fieldPoints );
  bsplineFilter->SetOrigin( this->m_Filter->m_Origin );
  bsplineFilter->SetSpacing( this->m_Filter->m_Spacing );
  bsplineFilter->SetSize( this->m_Filter->m_Size );
  bsplineFilter->SetDirection( this->m_Filter->m_Direction );
  bsplineFilter->SetNumberOfLevels( nlevels );
  bsplineFilter->SetSplineOrder( this->m_Filter->m_SplineOrder );
  bsplineFilter->SetNumberOfControlPoints( ncps );
  bsplineFilter->SetGenerateOutputImage( false );
  bsplineFilter->SetPointWeights( weights.GetPointer() );
  bsplineFilter->Update();

  derivative.SetSize( this->GetNumberOfParameters() );

  index = 0;
  ImageRegionIterator<typename TFilter::ControlPointLatticeType> ItP(
    bsplineFilter->GetPhiLattice(),
    bsplineFilter->GetPhiLattice()->GetLargestPossibleRegion() );
  for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
    {
    typename TFilter::ControlPointLatticeType::PixelType vector = ItP.Get();
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      derivative[index++] = vector[d];
      }
    }
}

template<class TFilter>
unsigned int
DMFFDLabeledPointSetRegistrationCostFunction<TFilter>
::GetNumberOfParameters() const
{
  unsigned numberOfParameters = Dimension * this->m_Filter->
    m_TotalDeformationFieldControlPointLattice->GetLargestPossibleRegion().
    GetNumberOfPixels();
  return numberOfParameters;
}

/**
 * DMFFDLabeledPointSetRegistrationFilter function definitions
 */
template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
DMFFDLabeledPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::DMFFDLabeledPointSetRegistrationFilter()
{

//  if( FixedPointSetType::Dimension !=
//         MovingPointSetType::Dimension )
//    {
//    itkExceptionMacro( "Point set dimensions must be equal." );
//    }

  this->SetNumberOfRequiredInputs( 2 );

  // Default values
  this->m_SplineOrder = 3;
  this->m_Size.Fill( 0 );
  this->m_Direction.SetIdentity();

  this->m_Optimizer = OptimizerType::New();
  this->m_MaximumNumberOfIterations.SetSize( 3 );
  this->m_MaximumNumberOfIterations.Fill( 10 );
  this->m_RelaxationFactor.Fill( 0.5 );

  this->m_UseInputAsSamples = true;
  this->m_NumberOfFixedSamples = 1000;
  this->m_NumberOfMovingSamples = 1000;
  this->m_UseAnisotropicCovariances = false;
  this->m_KernelSigma = 1.0;
  this->m_CovarianceKNeighborhood = 4;
  this->m_UseRegularizationTerm = true;
  this->m_AnnealingRate = 0.93;
  this->m_PointSetSigma = 1.0;
  this->m_Alpha = 1.0;
  this->m_EvaluationKNeighborhood = 50;

  this->m_Verbose = true;

  this->m_InitialDeformationFieldControlPointLattice
    = ControlPointLatticeType::New();
  this->m_InitialDeformationFieldControlPointLattice = NULL;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
DMFFDLabeledPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::~DMFFDLabeledPointSetRegistrationFilter()
{
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::GenerateData()
{
  for( this->m_CurrentLevel = 0; this->m_CurrentLevel
    < this->m_MaximumNumberOfIterations.Size(); this->m_CurrentLevel++ )
    {
    this->Initialize();

    if( this->m_Verbose )
      {
      std::cout << "Current level = " << this->m_CurrentLevel+1 << " ("
        << this->m_MaximumNumberOfIterations.Size() << " total levels).  "
        << "Control point grid size = "
        << this->m_TotalDeformationFieldControlPointLattice
          ->GetLargestPossibleRegion().GetSize() << ". " << std::endl;
      }

    itkDebugMacro( "Current level = " << this->m_CurrentLevel+1
      << " (" << this->m_MaximumNumberOfIterations.Size() << " total levels).  "
      << "Control point grid size = "
      << this->m_TotalDeformationFieldControlPointLattice
        ->GetLargestPossibleRegion().GetSize() << ". " );

    if( this->m_MaximumNumberOfIterations[this->m_CurrentLevel] > 0 )
      {
      typedef DMFFDLabeledPointSetRegistrationCostFunction<Self>
        RegistrationCostFunctionType;
      typename RegistrationCostFunctionType::Pointer costFunction =
        RegistrationCostFunctionType::New();
      costFunction->SetFilter( this );

      typename OptimizerType::ParametersType initialParameters;
      initialParameters.SetSize( costFunction->GetNumberOfParameters() );

      unsigned long index = 0;
      ImageRegionIterator<ControlPointLatticeType> It(
        this->m_TotalDeformationFieldControlPointLattice,
        this->m_TotalDeformationFieldControlPointLattice->
        GetLargestPossibleRegion() );
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        VectorType vector = It.Get();
        for( unsigned int d = 0; d < Dimension; d++ )
          {
          initialParameters[index++] = vector[d];
          }
        }

      RealType spacingNorm = 1.0;
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        spacingNorm *= vnl_math_sqr( this->m_Spacing[d] );
        }
      spacingNorm = vcl_sqrt( spacingNorm );

      typename OptimizerType::ScalesType scales;
      scales.SetSize( costFunction->GetNumberOfParameters() );
      scales.Fill( 1.0 );

      this->m_Optimizer->SetMinimize( true );
      this->m_Optimizer->SetCostFunction( costFunction );
      this->m_Optimizer->SetInitialPosition( initialParameters );
      this->m_Optimizer->SetScales( scales );

      if( this->m_CurrentLevel == 0 )
        {
        RealType stepLength = vcl_sqrt( this->m_LineSearchMaximumStepSize *
          spacingNorm *
          static_cast<RealType>( costFunction->GetNumberOfParameters() ) /
          static_cast<RealType>( Dimension ) );
        this->m_Optimizer->SetMaximumStepLength( stepLength );
        this->m_Optimizer->SetMinimumStepLength( 0.1 * stepLength );
        }
      else
        {
        this->m_Optimizer->SetMaximumStepLength(
          this->m_Optimizer->GetCurrentStepLength() );
        this->m_Optimizer->SetMinimumStepLength(
          this->m_Optimizer->GetMinimumStepLength() * 0.1 );
        }

      this->m_Optimizer->SetNumberOfIterations(
        this->m_MaximumNumberOfIterations[this->m_CurrentLevel] );
      this->m_Optimizer->SetRelaxationFactor(
        this->m_RelaxationFactor[this->m_CurrentLevel] );

      this->m_Optimizer->StartOptimization();

      typename OptimizerType::ParametersType
        finalParameters = this->m_Optimizer->GetCurrentPosition();

      index = 0;
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        VectorType vector = It.Get();
        for( unsigned int d = 0; d < Dimension; d++ )
          {
          vector[d] = finalParameters[index++];
          }
        It.Set( vector );
        }
      }
    }

  // Warp the input points

  typename WarpedPointSetType::Pointer warpedPoints
    = WarpedPointSetType::New();
  warpedPoints->Initialize();

  typename ControlPointFilterType::Pointer bspliner
    = ControlPointFilterType::New();
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPointLattice );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  bspliner->SetDirection( this->m_Direction );

  unsigned long count = 0;
  for( unsigned int j = 0; j < this->GetInput( 1 )->GetNumberOfPoints(); j++ )
    {
    MovingPointType inputPoint;
    this->GetInput( 1 )->GetPoint( j, &inputPoint );

    typename ControlPointFilterType::PointType point;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      point[d] = inputPoint[d];
      }
    VectorType vector;
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      continue;
      }
    point += vector;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      inputPoint[d] = point[d];
      }
    typename MovingPointSetType::PixelType inputLabel = 1;
    this->GetInput( 1 )->GetPointData( j, &inputLabel );

    warpedPoints->SetPoint( j, inputPoint );
    warpedPoints->SetPointData( count, inputLabel );
    count++;
    }
  this->SetNthOutput( 0, warpedPoints );
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::Initialize()
{
  if( this->m_CurrentLevel == 0 )
    {
    if( this->m_Verbose )
      {
      typename CommandIterationUpdate::Pointer observer =
        CommandIterationUpdate::New();
      this->m_Optimizer->AddObserver( itk::IterationEvent(), observer );
      }

    if( this->m_InitialDeformationFieldControlPointLattice )
      {
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        this->m_InitialMeshResolution[d] =
          this->m_InitialDeformationFieldControlPointLattice->
            GetLargestPossibleRegion().GetSize()[d] - this->m_SplineOrder;

        if( this->m_InitialMeshResolution[d] < 1 )
          {
          itkExceptionMacro( "Invalid size for initial deformation field "
            << "control point lattice." );
          }
        }
      }

    typename ControlPointLatticeType::RegionType::SizeType size;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      size[d] = this->m_InitialMeshResolution[d] + this->m_SplineOrder;
      }

    VectorType V;
    V.Fill( 0 );
    if( this->m_InitialDeformationFieldControlPointLattice )
      {
      typedef ImageDuplicator<ControlPointLatticeType> DuplicatorType;
      typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
      duplicator->SetInputImage( this->m_InitialDeformationFieldControlPointLattice );
      duplicator->Update();
      this->m_TotalDeformationFieldControlPointLattice = duplicator->GetOutput();
      }
    else
      {
      this->m_TotalDeformationFieldControlPointLattice
        = ControlPointLatticeType::New();
      this->m_TotalDeformationFieldControlPointLattice->SetRegions( size );
      this->m_TotalDeformationFieldControlPointLattice->Allocate();
      this->m_TotalDeformationFieldControlPointLattice->FillBuffer( V );
      }

    if( this->m_Size[0] == 0 )
      {

      /**
       * Define the transformation domain based on the bounding box of the
       * two point-sets.
       */
      FixedPointType minPoint;
      minPoint.Fill( NumericTraits<RealType>::max() );
      FixedPointType maxPoint;
      maxPoint.Fill( NumericTraits<RealType>::min() );
      for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
        {
        for( unsigned int d = 0; d < Dimension; d++ )
          {
          minPoint[d] = vnl_math_min( minPoint[d],
            this->GetInput( i )->GetBoundingBox()->GetMinimum()[d] );
          maxPoint[d] = vnl_math_max( maxPoint[d],
            this->GetInput( i )->GetBoundingBox()->GetMaximum()[d] );
          }
        }

      this->m_Size.Fill( 100 );

      for( unsigned int d = 0; d < Dimension; d++ )
        {
        this->m_Origin[d] = minPoint[d];
        this->m_Spacing[d] = ( maxPoint[d] - minPoint[d] )
          / static_cast<RealType>( this->m_Size[d] - 1 );
        }
      }
    }
  else
    {
    typename ControlPointFilterType::Pointer bspliner
      = ControlPointFilterType::New();

    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetInput( this->m_TotalDeformationFieldControlPointLattice );

    typename BSplineFilterType::ArrayType nlevels;
    nlevels.Fill( 2 );
    this->m_TotalDeformationFieldControlPointLattice
      = bspliner->RefineControlPointLattice( nlevels );
    }
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Transformation variables" << std::endl;
  os << indent << "   " << "Spline order: "
    << this->m_SplineOrder << std::endl;
  os << indent << "  " << "Initial mesh resolution: "
    << this->m_InitialMeshResolution << std::endl;
  os << indent << "  " << "Origin: "
    << this->m_Origin << std::endl;
  os << indent << "  " << "Spacing: "
    << this->m_Spacing << std::endl;
  os << indent << "  " << "Size: "
    << this->m_Size << std::endl;

  os << indent << "Jensen-Havrda-Charvat-Tsallis similarity variables"
    << std::endl;
  os << indent << "Use regularization term: "
     << this->m_UseRegularizationTerm << std::endl;
  os << indent << "Alpha: "
     << this->m_Alpha << std::endl;

  os << indent << "Point-set sigma: "
     << this->m_PointSetSigma << std::endl;
  if( !this->m_UseInputAsSamples )
    {
    os << indent << "Number of fixed samples: "
       << this->m_NumberOfFixedSamples << std::endl;
    os << indent << "Number of moving samples: "
       << this->m_NumberOfMovingSamples << std::endl;
    }
  else
    {
    os << indent << "Use input points as samples." << std::endl;
    }

  if( this->m_UseAnisotropicCovariances )
    {
    os << indent << "Kernel sigma: "
       << this->m_KernelSigma << std::endl;
    os << indent << "Covariance k-neighborhood: "
       << this->m_CovarianceKNeighborhood << std::endl;
    }
  else
    {
    os << indent << "Isotropic covariances are used." << std::endl;
    }
}

}
#endif
