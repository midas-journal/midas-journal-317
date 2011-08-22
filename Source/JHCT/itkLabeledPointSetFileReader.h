/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabeledPointSetFileReader.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:21:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabeledPointSetFileReader_h
#define __itkLabeledPointSetFileReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"

#include "itkImage.h"

#include <vector>

namespace itk
{

/** \class LabeledPointSetFileReader
 * \brief
 * Reads a file and creates an itkMesh.
 *
 */
template <class TOutputMesh>
class LabeledPointSetFileReader 
: public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef LabeledPointSetFileReader               Self;
  typedef MeshSource<TOutputMesh>                 Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TOutputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( LabeledPointSetFileReader, MeshSource );

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                             OutputMeshType;
  typedef typename OutputMeshType::MeshTraits     MeshTraits;
  typedef typename OutputMeshType::Superclass     PointSetType;
  typedef typename OutputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType          PixelType;

  typedef Image<PixelType, 
    itkGetStaticConstMacro( Dimension )>          LabeledPointSetImageType;

  typedef std::vector<PixelType>                  LabelSetType;

  /** Set/Get the name of the file to be read. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );   
  
  itkSetMacro( ExtractBoundaryPoints, bool );
  itkGetMacro( ExtractBoundaryPoints, bool );
  itkBooleanMacro( ExtractBoundaryPoints );

  /**
   * Percentage of points selected randomnly
   */
  itkSetClampMacro( RandomPercentage, double, 0.0, 1.0 );
  itkGetConstMacro( RandomPercentage, double );

  LabelSetType* GetLabelSet() { return &this->m_LabelSet; }  
  unsigned int GetNumberOfLabels() { return this->m_LabelSet.size(); }  

protected:
  LabeledPointSetFileReader();
  ~LabeledPointSetFileReader() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** Reads the file */
  void GenerateData();

  bool                                            m_ExtractBoundaryPoints;

  std::string                                     m_FileName;  
  double                                          m_RandomPercentage;
  LabelSetType                                    m_LabelSet;
  

private:
  LabeledPointSetFileReader( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented
  
  void ReadLandmarksFromImageFile();
  void ReadLandmarksFromAvantsFile();
  void ReadLandmarksFromVTKFile();

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabeledPointSetFileReader.txx"
#endif

#endif //_itkLabeledPointSetFileReader_h
