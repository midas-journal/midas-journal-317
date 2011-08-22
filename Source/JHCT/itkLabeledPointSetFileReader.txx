/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabeledPointSetFileReader.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:21:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabeledPointSetFileReader_txx
#define __itkLabeledPointSetFileReader_txx

#include "itkLabeledPointSetFileReader.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include <fstream>
#include <stdio.h>
#include <string>

namespace itk
{

//
// Constructor
//
template<class TOutputMesh>
LabeledPointSetFileReader<TOutputMesh>
::LabeledPointSetFileReader()
{
  this->m_RandomPercentage = 1.0;
  this->m_ExtractBoundaryPoints = false;
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}

template<class TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No input FileName" );
    return;
    }

  //
  // Read input file
  //
  std::ifstream inputFile( m_FileName.c_str() );

  if( !inputFile.is_open() )
    {
    itkExceptionMacro("Unable to open file\n"
        "inputFilename= " << m_FileName );
    return;
    }
  else
    {
    inputFile.close();
    }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "txt" )
    {
    this->ReadLandmarksFromAvantsFile();
    }
  else if( extension == "vtk" )
    {
    this->ReadLandmarksFromVTKFile();
    }
  else // try reading the file as an image
    {
    this->ReadLandmarksFromImageFile();
    }

  if( this->m_RandomPercentage < 1.0 )
    {
    typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    typename GeneratorType::Pointer generator = GeneratorType::New();
    generator->SetSeed();

    typename OutputMeshType::Pointer output = OutputMeshType::New();
    output->Initialize();
    
    if( this->GetOutput()->GetNumberOfPoints() > 0 )
      {
      typename OutputMeshType::PointsContainerIterator It =
        this->GetOutput()->GetPoints()->Begin();
      typename OutputMeshType::PointDataContainerIterator ItD =
        this->GetOutput()->GetPointData()->Begin();
  
      unsigned long count = 0;
      while( It != this->GetOutput()->GetPoints()->End() )
        {
        if( generator->GetVariateWithClosedRange() <= this->m_RandomPercentage )
          {
          output->SetPoint( count, It.Value() );
          output->SetPointData( count, ItD.Value() );
          count++;
          }
        ++It;
        ++ItD;
        }
      }  
    this->GraftOutput( output );
    }
  this->m_LabelSet.clear();

  if( this->GetOutput()->GetNumberOfPoints() > 0 )
    {
    typename OutputMeshType::PointDataContainerIterator ItD =
      this->GetOutput()->GetPointData()->Begin();
    while( ItD != this->GetOutput()->GetPointData()->End() )
      {
      if( find( this->m_LabelSet.begin(), this->m_LabelSet.end(), ItD.Value() )
             == this->m_LabelSet.end() )
        {
        this->m_LabelSet.push_back( ItD.Value() );
        }
      ++ItD;
      }
    }
}

template<class TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>
::ReadLandmarksFromAvantsFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile( m_FileName.c_str() );

  unsigned long count = 0;
  while( !inputFile.eof() )
    {
    PointType point;
    PixelType label;

    if( Dimension == 2 )
      {
      float trash; 
      inputFile >> point >> trash >> label;
      }
    else // Dimension == 3
      {
      inputFile >> point >> label;
      }  
    if( ( point.GetVectorFromOrigin() ).GetSquaredNorm() > 0.0
         || label != 0 )
      {
      outputMesh->SetPointData( count, label );
      outputMesh->SetPoint( count, point );
      count++;
      }
    }
  inputFile.close();
}

template<class TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>
::ReadLandmarksFromVTKFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile( m_FileName.c_str() );

  std::string line;

  while( !inputFile.eof() )
    {
    std::getline( inputFile, line );

    if( line.find( "POINTS" ) != std::string::npos )
      {
      break;
      }
    }

  itkDebugMacro( "POINTS line" << line );

  std::string pointLine( line, strlen( "POINTS " ), line.length() );
  itkDebugMacro( "pointLine " << pointLine );

  int numberOfPoints = -1;

  if( sscanf( pointLine.c_str(),"%d",&numberOfPoints ) != 1 )
    {
    itkExceptionMacro( "ERROR: Failed to read numberOfPoints\n"
        "       pointLine = " << pointLine );
    return;
    }

  itkDebugMacro( "numberOfPoints = " << numberOfPoints );

  if( numberOfPoints < 1 )
    {
    itkExceptionMacro( "numberOfPoints < 1"
        << "       numberOfPoints = " << numberOfPoints );
    return;
    }

  outputMesh->GetPoints()->Reserve( numberOfPoints );

  //
  // Load the point coordinates into the itk::Mesh
  //
  PointType point;

  for( long i = 0; i < numberOfPoints; i++ )
    {
    if( Dimension == 2 ) 
      {
      float trash; 
      inputFile >> point >> trash; 
      }
    else  // Dimension = 3
      { 
      inputFile >> point;
      }  
    outputMesh->SetPoint( i, point );
    }

  //
  // Find the labels associated with each pixel
  //
  while( !inputFile.eof() )
    {
    std::getline( inputFile, line );

    if( line.find( "SCALARS" ) != std::string::npos )
      {
      break;
      }
    }

  if( inputFile.eof() )
    {
    inputFile.close();
    return;
    }

  std::string::size_type pos = line.rfind( " " );
  std::string numberOfComponents( line, pos+1, line.length()-1 );

  if( numberOfComponents != "1" )
    {
    itkExceptionMacro( "Only single label components are readable" );
    }
  std::getline( inputFile, line );

  PixelType label;
  for( long i = 0; i < numberOfPoints; i++ )
    {
    inputFile >> label;
    outputMesh->SetPointData( i, label );
    }

  inputFile.close();
}

template<class TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>
::ReadLandmarksFromImageFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  typedef ImageFileReader<LabeledPointSetImageType> ImageReaderType;
  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( this->m_FileName.c_str() );
  imageReader->Update();

  if( !this->m_ExtractBoundaryPoints )
    {
    ImageRegionIteratorWithIndex<LabeledPointSetImageType> It( imageReader->GetOutput(),
      imageReader->GetOutput()->GetLargestPossibleRegion() );
  
    unsigned long count = 0;
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( label != NumericTraits<PixelType>::Zero )
        {
        typename LabeledPointSetImageType::PointType imagePoint;
        imageReader->GetOutput()->TransformIndexToPhysicalPoint(
          It.GetIndex(), imagePoint );
  
        PointType point;
        point.CastFrom( imagePoint );
        outputMesh->SetPoint( count, point );
        outputMesh->SetPointData( count, label );
        count++;
        }
      }
    }
  else
    {    
    this->m_LabelSet.clear();
  
    ImageRegionIterator<LabeledPointSetImageType> It ( imageReader->GetOutput(),
      imageReader->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType label = It.Get();
      if( It.Get() > 0 )
        {
        if( find( this->m_LabelSet.begin(), this->m_LabelSet.end(), label )
               == this->m_LabelSet.end() )
          {
          this->m_LabelSet.push_back( label );
          }
        }  
      }

    unsigned long count = 0;

    typename LabelSetType::const_iterator it;
    for( it = this->m_LabelSet.begin(); it != this->m_LabelSet.end(); ++it )
      {
      typedef Image< unsigned char, Dimension> BinaryImageType;
      typedef BinaryThresholdImageFilter<LabeledPointSetImageType, 
        BinaryImageType> BinaryFilterType;
      typename BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();
    
      binaryFilter->SetLowerThreshold( *it );
      binaryFilter->SetUpperThreshold( *it );
      binaryFilter->SetInsideValue( 0 );
      binaryFilter->SetOutsideValue( 1 );
      binaryFilter->SetInput( imageReader->GetOutput() );
      binaryFilter->Update();
    
      typedef BinaryBallStructuringElement<unsigned char,
        Dimension> StructuringElementType;
    
      typedef BinaryErodeImageFilter<BinaryImageType, BinaryImageType,
                              StructuringElementType> ErodeType;
      typename ErodeType::Pointer erode = ErodeType::New();
    
      StructuringElementType structuringElement;
      structuringElement.SetRadius( 1 );
      structuringElement.CreateStructuringElement();
      erode->SetKernel( structuringElement );
      erode->SetForegroundValue( 1 );
      erode->SetBackgroundValue( 2 );
      erode->SetInput( binaryFilter->GetOutput() );
      erode->Update();
      
      ImageRegionIterator<BinaryImageType> ItE( erode->GetOutput(),
        erode->GetOutput()->GetLargestPossibleRegion() );

      for( ItE.GoToBegin(); !ItE.IsAtEnd(); ++ItE )
        {
        if( ItE.Get() == 2 )
          {
          typename LabeledPointSetImageType::PointType imagePoint;
          imageReader->GetOutput()->TransformIndexToPhysicalPoint(
            ItE.GetIndex(), imagePoint );
    
          PointType point;
          point.CastFrom( imagePoint );
          outputMesh->SetPoint( count, point );
          outputMesh->SetPointData( count, *it );
          count++;
          }  
        }
      }
    }  
}

template<class TOutputMesh>
void
LabeledPointSetFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << this->m_FileName << std::endl;
  os << indent << "RandomPercentage: "
     << this->m_RandomPercentage << std::endl;

}

} //end of namespace itk


#endif
