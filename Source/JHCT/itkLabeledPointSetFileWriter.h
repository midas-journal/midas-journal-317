/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabeledPointSetFileWriter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/05 16:09:11 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabeledPointSetFileWriter_h
#define __itkLabeledPointSetFileWriter_h

#include "itkMesh.h"
#include "itkTriangleCell.h"

namespace itk
{
/** \class LabeledPointSetFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TInputMesh>
class LabeledPointSetFileWriter : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef LabeledPointSetFileWriter                      Self;
  typedef Object                                  Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void );
  void Write( void );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TInputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( LabeledPointSetFileWriter, Object );

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputMesh                             InputMeshType;
  typedef typename TInputMesh::Pointer           InputMeshPointer;
  typedef typename InputMeshType::MeshTraits     MeshTraits;
  typedef typename InputMeshType::Superclass     PointSetType;
  typedef typename InputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType         PixelType;
  typedef Image<PixelType,
    itkGetStaticConstMacro( Dimension )>         LabeledPointSetImageType;
  typedef typename
    LabeledPointSetImageType::SizeType           ImageSizeType;
  typedef typename
    LabeledPointSetImageType::PointType          ImageOriginType;
  typedef typename
    LabeledPointSetImageType::SpacingType        ImageSpacingType;
  typedef typename
    LabeledPointSetImageType::DirectionType      ImageDirectionType;


  /** Set the Input */
  void SetInput( InputMeshType * input );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  /** Specify image attributes if output is an image. */
  itkSetMacro( ImageSize, ImageSizeType );
  itkGetConstMacro( ImageSize, ImageSizeType );

  itkSetMacro( ImageOrigin, ImageOriginType );
  itkGetConstMacro( ImageOrigin, ImageOriginType );

  itkSetMacro( ImageSpacing, ImageSpacingType );
  itkGetConstMacro( ImageSpacing, ImageSpacingType );

  itkSetMacro( ImageDirection, ImageDirectionType );
  itkGetConstMacro( ImageDirection, ImageDirectionType );

protected:
  LabeledPointSetFileWriter();
  virtual ~LabeledPointSetFileWriter();

  virtual void GenerateData();

  std::string                         m_FileName;
  InputMeshPointer                    m_Input;

  /**
   * If output is an image type, the attributes must be specified.
   */
  ImageSizeType                       m_ImageSize;
  ImageSpacingType                    m_ImageSpacing;
  ImageOriginType                     m_ImageOrigin;
  ImageDirectionType                  m_ImageDirection;

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  LabeledPointSetFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void WriteLandmarksToAvantsFile();
  void WriteLandmarksToVTKFile();
  void WriteLandmarksToImageFile();

};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabeledPointSetFileWriter.txx"
#endif

#endif
