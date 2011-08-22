/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabeledPointSetFileWriter.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/05 16:09:11 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabeledPointSetFileWriter_txx
#define __itkLabeledPointSetFileWriter_txx

#include "itkLabeledPointSetFileWriter.h"
#include "itkImageFileWriter.h"

#include <fstream>

namespace itk
{

//
// Constructor
//
template<class TInputMesh>
LabeledPointSetFileWriter<TInputMesh>
::LabeledPointSetFileWriter()
{
  this->m_Input = NULL;
  this->m_FileName = "";

  this->m_ImageSize.Fill( 0 );
}

//
// Destructor
//
template<class TInputMesh>
LabeledPointSetFileWriter<TInputMesh>
::~LabeledPointSetFileWriter()
{
}

//
// Set the input mesh
//
template<class TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>
::SetInput(InputMeshType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void LabeledPointSetFileWriter<TInputMesh>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void LabeledPointSetFileWriter<TInputMesh>
::Write()
{
  this->GenerateData();
}

template<class TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No FileName" );
    return;
    }

  if( this->m_ImageSize[0] == 0 )
    {
    this->m_ImageSize.Fill( 100 );
    this->m_ImageOrigin.CastFrom( 
      this->m_Input->GetBoundingBox()->GetMinimum() );
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_ImageSpacing[d] = ( 
        this->m_Input->GetBoundingBox()->GetMaximum()[d] -
        this->m_Input->GetBoundingBox()->GetMinimum()[d] ) 
        / static_cast<double>( this->m_ImageSize[d] + 1 );
      }
    this->m_ImageDirection.SetIdentity();
    }

  //
  // Read output file
  //
  std::ofstream outputFile( m_FileName.c_str() );

  if( !outputFile.is_open() )
    {
    itkExceptionMacro( "Unable to open file\n"
        "outputFilename= " << m_FileName );
    return;
    }
  else
    {
    outputFile.close();
    }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "txt" )
    {
    this->WriteLandmarksToAvantsFile();
    }
  else if( extension == "vtk" )
    {
    this->WriteLandmarksToVTKFile();
    }
  else
    {
    try
      {
      this->WriteLandmarksToImageFile();
      }
    catch(...)
      {
      itkExceptionMacro( "Unknown extension: " << extension );
      }
    }

}

template<class TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>
::WriteLandmarksToVTKFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str() );

  outputFile << "# vtk DataFile Version 2.0" << std::endl;
  outputFile << "File written by itkLabeledPointSetFileWriter" << std::endl;
  outputFile << "ASCII" << std::endl;
  outputFile << "DATASET POLYDATA" << std::endl;

  // POINTS go first

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();
  outputFile << "POINTS " << numberOfPoints << " float" << std::endl;

  typename InputMeshType::PointsContainerIterator pointIterator
    = this->m_Input->GetPoints()->Begin();
  typename InputMeshType::PointsContainerIterator pointEnd
    = this->m_Input->GetPoints()->End();
  while( pointIterator != pointEnd )
    {
    PointType point = pointIterator.Value();
    outputFile << point[0] << " " << point[1];
    if( Dimension == 2 )
      {
      outputFile << " 0 " << std::endl;
      }
    else if( Dimension == 3 )
      {
      outputFile << " " << point[2] << " " << std::endl;
      }
    pointIterator++;
    }

  outputFile << std::endl;
  outputFile << "POINT_DATA " << numberOfPoints << std::endl;
  outputFile << "SCALARS pointLabels long 1" << std::endl;
  outputFile << "LOOKUP_TABLE default" << std::endl;

  typename InputMeshType::PointDataContainerIterator pointDataIterator
    = this->m_Input->GetPointData()->Begin();
  typename InputMeshType::PointDataContainerIterator pointDataEnd
    = this->m_Input->GetPointData()->End();

  while( pointDataIterator != pointDataEnd )
    {
    outputFile << pointDataIterator.Value() << " ";
    pointDataIterator++;
    }
  outputFile << std::endl;

  outputFile.close();
}

template<class TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>
::WriteLandmarksToAvantsFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str() );

  outputFile << "0 0 0 0" << std::endl;

  if( this->m_Input->GetNumberOfPoints() > 0 )
    {
    typename InputMeshType::PointsContainerIterator pointIterator
      = this->m_Input->GetPoints()->Begin();
    typename InputMeshType::PointsContainerIterator pointEnd
      = this->m_Input->GetPoints()->End();

    typename InputMeshType::PointDataContainerIterator pointDataIterator
      = this->m_Input->GetPointData()->Begin();

    while( pointIterator != pointEnd )
      {
      PointType point = pointIterator.Value();
      outputFile << point[0] << " " << point[1];
      if( Dimension == 2 )
        {
        outputFile << " 0 ";
        }
      else if( Dimension == 3 )
        {
        outputFile << " " << point[2] << " ";
        }
      outputFile << pointDataIterator.Value() << std::endl;
      pointIterator++;
      pointDataIterator++;
      }
    }

  outputFile << "0 0 0 0" << std::endl;

  outputFile.close();
}

template<class TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>
::WriteLandmarksToImageFile()
{
  typename LabeledPointSetImageType::Pointer outputImage
    = LabeledPointSetImageType::New();
  outputImage->SetDirection( this->m_ImageDirection );
  outputImage->SetRegions( this->m_ImageSize );
  outputImage->SetOrigin( this->m_ImageOrigin );
  outputImage->SetSpacing( this->m_ImageSpacing );
  outputImage->Allocate();
  outputImage->FillBuffer( NumericTraits<PixelType>::Zero );

  if( this->m_Input->GetNumberOfPoints() > 0 )
    {
    typename InputMeshType::PointsContainerIterator pointIterator
      = this->m_Input->GetPoints()->Begin();
    typename InputMeshType::PointsContainerIterator pointEnd
      = this->m_Input->GetPoints()->End();

    typename InputMeshType::PointDataContainerIterator pointDataIterator
      = this->m_Input->GetPointData()->Begin();

    while( pointIterator != pointEnd )
      {
      PointType point = pointIterator.Value();
      PixelType label = pointDataIterator.Value();

      typename LabeledPointSetImageType::IndexType index;
      typename LabeledPointSetImageType::PointType ipoint;
      ipoint.CastFrom( point );
      if( outputImage->TransformPhysicalPointToIndex( ipoint, index ) )
        {
        outputImage->SetPixel( index, label );
        }
      pointIterator++;
      pointDataIterator++;
      }
    }

  typedef ImageFileWriter<LabeledPointSetImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( this->m_FileName.c_str() );
  writer->SetInput( outputImage );
  writer->Update();
}

template<class TInputMesh>
void
LabeledPointSetFileWriter<TInputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
