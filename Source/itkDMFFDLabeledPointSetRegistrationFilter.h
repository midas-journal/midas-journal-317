/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDMFFDLabeledPointSetRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDMFFDLabeledPointSetRegistrationFilter_h
#define __itkDMFFDLabeledPointSetRegistrationFilter_h

#include "itkPointSetToPointSetFilter.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkPointSet.h"
#include "itkSingleValuedCostFunction.h"
#include "itkVector.h"

#include <map>

/**
 * \class DMFFDLabeledPointSetRegistrationFilter
 * \brief Point set registration filter using the Jensen Havrda Charvat Tsallis
 * entropy and a directly manipulated free-form deformation transformation
 * model.
 *
 * \par Information-theoretic labeled point-set registration filter which
 * uses a directly manipulated free-form deformation model.  Input consists
 * of two labeled point-sets (fixed and moving).  Output is the warped
 * point-set.  The total deformation field, described by control points,
 * can subsequently be used to generate the resulting sampled deformation
 * field.
 *
 * \par REFERENCE
 * NJ Tustison, SP Awate, G Song, TS Cook, and JC Gee, "A new information-
 * theoretic measure to control the robustness-sensitivity trade-off for
 * DMFFD point-set registration." Proceedings of the 21st Biennial
 * International Conference on Information Processing in Medical Imaging
 * (IPMI), 2009
*/

namespace itk {

class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate   Self;
  itkNewMacro( Self );

  typedef RegularStepGradientDescentOptimizer OptimizerType;

protected:
  CommandIterationUpdate() {};

public:
  void Execute( itk::Object *caller, const itk::EventObject & event )
    {
    Execute( ( const itk::Object * ) caller, event );
    }

  void Execute( const itk::Object * object, const itk::EventObject & event )
    {
    const OptimizerType * optimizer =
      dynamic_cast<const OptimizerType *>( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << "   Iteration " << optimizer->GetCurrentIteration();
    std::cout << " (step length = "
      << optimizer->GetCurrentStepLength() << "): ";
    std::cout << optimizer->GetValue() << std::endl;
    }
};

/**
 * Class definition for N3BiasFieldScaleCostFunction
 */

template<class TFilter>
class ITK_EXPORT DMFFDLabeledPointSetRegistrationCostFunction
  : public SingleValuedCostFunction
{
public:
  typedef DMFFDLabeledPointSetRegistrationCostFunction   Self;
  typedef SingleValuedCostFunction                       Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( DMFFDLabeledPointSetRegistrationCostFunction,
    SingleValuedCostFunction );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef typename TFilter::RealType     RealType;

  typedef Superclass::MeasureType        MeasureType;
  typedef Superclass::DerivativeType     DerivativeType;
  typedef Superclass::ParametersType     ParametersType;

  itkSetMacro( Filter, typename TFilter::Pointer );
  itkGetConstMacro( Filter, typename TFilter::Pointer );

  itkStaticConstMacro( Dimension, unsigned int,
    TFilter::FixedPointSetType::PointDimension );

  virtual MeasureType GetValue( const ParametersType & parameters ) const;
  virtual void GetDerivative( const ParametersType & parameters,
    DerivativeType & derivative ) const;
  virtual void GetValueAndDerivative( const ParametersType & parameters,
    MeasureType & value, DerivativeType & derivative ) const;
  virtual unsigned int GetNumberOfParameters() const;

protected:
  DMFFDLabeledPointSetRegistrationCostFunction();
  virtual ~DMFFDLabeledPointSetRegistrationCostFunction();

private:
  DMFFDLabeledPointSetRegistrationCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename TFilter::Pointer              m_Filter;
};

template<class TFixedPointSet, class TMovingPointSet = TFixedPointSet,
  class TWarpedPointSet = TFixedPointSet>
class ITK_EXPORT DMFFDLabeledPointSetRegistrationFilter
: public PointSetToPointSetFilter<TMovingPointSet, TWarpedPointSet>
{
public:
  typedef DMFFDLabeledPointSetRegistrationFilter         Self;
  typedef PointSetToPointSetFilter
    <TMovingPointSet, TWarpedPointSet>                   Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;


  template<class Self> friend class DMFFDLabeledPointSetRegistrationCostFunction;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( DMFFDLabeledPointSetRegistrationFilter,
    PointSetToPointSetFilter );

  /**
   * Point set typedefs
   */
  typedef TFixedPointSet                            FixedPointSetType;
  typedef typename FixedPointSetType::PointType     FixedPointType;
  typedef typename FixedPointSetType::PixelType     FixedPixelType;
  typedef TMovingPointSet                           MovingPointSetType;
  typedef typename MovingPointSetType::PointType    MovingPointType;
  typedef typename MovingPointSetType::PixelType    MovingPixelType;
  typedef TWarpedPointSet                           WarpedPointSetType;
  typedef typename WarpedPointSetType::PointType    WarpedPointType;
  typedef typename WarpedPointSetType::PixelType    WarpedPixelType;

  /**
   * Dimensionality of input and output data is assumed to be the same.
   */
  itkStaticConstMacro( Dimension, unsigned int,
    FixedPointSetType::PointDimension );

  typedef float                                     RealType;
  typedef Image<RealType,
    itkGetStaticConstMacro( Dimension )>            RealImageType;
  typedef FixedArray<unsigned int,
    itkGetStaticConstMacro( Dimension )>            ArrayType;
  typedef Array<unsigned int>                       ResizableUIntArrayType;
  typedef Array<RealType>                           ResizableRealArrayType;
  typedef std::map<WarpedPixelType, RealType>       LabelWeightsMapType;

  typedef RegularStepGradientDescentOptimizer      OptimizerType;

  /** Typedefs for B-spline filter */
  typedef Vector<RealType,
    itkGetStaticConstMacro( Dimension )>            VectorType;
  typedef Image<VectorType,
    itkGetStaticConstMacro( Dimension )>            DeformationFieldType;
  typedef PointSet<VectorType,
    itkGetStaticConstMacro( Dimension )>            BSplinePointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <BSplinePointSetType, DeformationFieldType>     BSplineFilterType;
  typedef typename
    BSplineFilterType::PointDataImageType           ControlPointLatticeType;
  typedef typename ControlPointLatticeType::Pointer ControlPointLatticePointer;
  typedef typename
    BSplineFilterType::WeightsContainerType         WeightsContainerType;
  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType> ControlPointFilterType;
  typedef typename
    ControlPointFilterType::OriginType              OriginType;
  typedef typename
    ControlPointFilterType::SpacingType             SpacingType;
  typedef typename ControlPointFilterType::SizeType SizeType;
  typedef typename
    ControlPointFilterType::DirectionType           DirectionType;

  /** Set/Get functions */

  /**
   * Variables governing the B-spline transformation model
   */

  /**
   * B-spline order
   */
  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  /**
   * B-spline mesh size at the lowest resolution.
   */
  itkSetMacro( InitialMeshResolution, ArrayType );
  itkGetConstMacro( InitialMeshResolution, ArrayType );

  /**
   * Weights for each label to emphasize different labels in the gradient.
   */
  void SetLabelWeights( LabelWeightsMapType labelWeights )
    {
    this->m_LabelWeights = labelWeights;
    this->Modified();
    }

  /**
   * The origin, spacing, and size define the transformation domain.
   */
  itkSetMacro( Origin, OriginType );
  itkGetConstMacro( Origin, OriginType );

  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );

  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );

  itkSetMacro( Direction, DirectionType );
  itkGetConstMacro( Direction, DirectionType );

  /**
   * Variables governing the gradient descent iterative optimization.
   */

  /**
   * Maximum number of iterations for each resolution level.
   */
  itkSetMacro( MaximumNumberOfIterations, ResizableUIntArrayType );
  itkGetConstMacro( MaximumNumberOfIterations, ResizableUIntArrayType );

  /**
   * Relaxation factor.
   */
  itkSetMacro( RelaxationFactor, ResizableRealArrayType );
  itkGetConstMacro( RelaxationFactor, ResizableRealArrayType );

  /**
   * Maximum step size for each line iteration.
   */
  itkSetMacro( LineSearchMaximumStepSize, RealType );
  itkGetConstMacro( LineSearchMaximumStepSize, RealType );

  /**
   * Variables governing the Jensen-Havrda-Charvat-Tsallis point-set
   * similarity measure.
   */

  /**
   * To calculate the gradient, one can either set the number of fixed
   * and moving samples to be generated randomly from the input point-sets or
   * one can use the input point-sets as samples.
   */
  itkSetMacro( NumberOfFixedSamples, unsigned long );
  itkGetConstMacro( NumberOfFixedSamples, unsigned long );

  itkSetMacro( NumberOfMovingSamples, unsigned long );
  itkGetConstMacro( NumberOfMovingSamples, unsigned long );

  itkSetMacro( UseInputAsSamples, bool );
  itkGetConstMacro( UseInputAsSamples, bool );
  itkBooleanMacro( UseInputAsSamples );

  /**
   * Anisotropic covariances are calculated from the local point-set
   * structure (using m_KernelSigma and m_Covariance KNeighborhood as a
   * weighting parameter) and the covariance generated from
   * m_PointSetSigma. and the current annealing value  Otherwise,
   * isotropic covariances are used employing m_PointSetSigma and the
   * current annealing value.
   */
  itkSetMacro( UseAnisotropicCovariances, bool );
  itkGetConstMacro( UseAnisotropicCovariances, bool );
  itkBooleanMacro( UseAnisotropicCovariances );

  itkSetMacro( KernelSigma, RealType );
  itkGetConstMacro( KernelSigma, RealType );

  itkSetMacro( CovarianceKNeighborhood, unsigned int );
  itkGetConstMacro( CovarianceKNeighborhood, unsigned int );

  itkSetMacro( PointSetSigma, RealType );
  itkGetConstMacro( PointSetSigma, RealType );

  /**
   * Describes the deterministic annealing schedule.
   */
  itkSetClampMacro( AnnealingRate, RealType, 0.00001, 1.0 );
  itkGetMacro( AnnealingRate, RealType );

  /**
   * Use the second term of the JHCT divergence.
   */
  itkSetMacro( UseRegularizationTerm, bool );
  itkGetConstMacro( UseRegularizationTerm, bool );
  itkBooleanMacro( UseRegularizationTerm );

  /**
   * The alpha parameter of the JHCT divergence controls the robustness
   * sensitivity trade-off with alpha = 1 equivalent to the ML solution and
   * alpha = 2 the robust L2 solution.
   */
  itkSetMacro( Alpha, RealType );
  itkGetConstMacro( Alpha, RealType );

  /**
   * Specifies the number of nearest neighbors (k-d tree) to be used in
   * calculating the metric and gradient values.
   */
  itkSetMacro( EvaluationKNeighborhood, unsigned int );
  itkGetConstMacro( EvaluationKNeighborhood, unsigned int );

  /**
   * Print to screen the course of optimization.
   */
  itkSetMacro( Verbose, bool );
  itkGetConstMacro( Verbose, bool );
  itkBooleanMacro( Verbose );

  /**
   * Allow one to specify and initial deformation field defined by a
   * control point lattic.e
   */
  itkSetObjectMacro( InitialDeformationFieldControlPointLattice,
    ControlPointLatticeType );
  itkGetObjectMacro( InitialDeformationFieldControlPointLattice,
    ControlPointLatticeType );

  /**
   * Get the resulting transformation model defined by a grid of control point
   * grids.
   */
  itkGetObjectMacro( TotalDeformationFieldControlPointLattice,
    ControlPointLatticeType );

protected :
  /** de/constructor */
  DMFFDLabeledPointSetRegistrationFilter();
  ~DMFFDLabeledPointSetRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :
  //purposely not implemented
  DMFFDLabeledPointSetRegistrationFilter( const Self& );
  void operator=( const Self& ); //purposely not implemented

  void Initialize();

private :

  /**
   * Variables governing the B-spline transformation model
   */
  unsigned int                          m_SplineOrder;
  ArrayType                             m_InitialMeshResolution;
  LabelWeightsMapType                   m_LabelWeights;
  OriginType                            m_Origin;
  SpacingType                           m_Spacing;
  SizeType                              m_Size;
  DirectionType                         m_Direction;

  /**
   * Variables governing the  optimization
   */
  typename OptimizerType::Pointer       m_Optimizer;
  ResizableUIntArrayType                m_MaximumNumberOfIterations;
  RealType                              m_LineSearchMaximumStepSize;
  ResizableRealArrayType                m_RelaxationFactor;
  unsigned int                          m_CurrentLevel;
  int                                   m_CurrentIteration;

  /**
   * Variables governing the Jensen-Havrda-Charvat-Tsallis point-set
   * similarity measure
   */
  bool                                  m_UseInputAsSamples;
  unsigned long                         m_NumberOfFixedSamples;
  unsigned long                         m_NumberOfMovingSamples;
  bool                                  m_UseAnisotropicCovariances;
  RealType                              m_KernelSigma;
  unsigned int                          m_CovarianceKNeighborhood;
  RealType                              m_PointSetSigma;
  RealType                              m_AnnealingRate;
  bool                                  m_UseRegularizationTerm;
  RealType                              m_Alpha;
  unsigned int                          m_EvaluationKNeighborhood;

  /**
   * Control point lattices
   */
  ControlPointLatticePointer            m_InitialDeformationFieldControlPointLattice;
  ControlPointLatticePointer            m_TotalDeformationFieldControlPointLattice;

  /**
   * Other variables.
   */
  bool                                  m_Verbose;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDMFFDLabeledPointSetRegistrationFilter.txx"
#endif

#endif
