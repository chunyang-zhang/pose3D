/*
Cameras.h

Motion class  

1. read an AMC file and store it in a sequence of state vector 
2. write an AMC file
3. export to a mrdplot format for plotting the trajectories

You can add more motion data processing functions in this class. 

Revision 1 - Steve Lin (CMU), Jan 14, 2002
Revision 2 - Alla Safonova and Kiran Bhat (CMU), Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao (USC), Feb, 2012

*/
#ifndef _CAMERAS_H_
#define _CAMERAS_H_

#include <FL/glu.h>
#include"types.h"
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include"skeleton.h"



struct PinholeCamera{
  int poseId;
  double w;
  double h;
  Eigen::Matrix3d K;
  double sz;
  std::vector<Eigen::Matrix4d> pos;    
};

class Cameras 
{
  //function members
public:    
  Cameras();
  ~Cameras();  
  int loadNewCameraFromFile(char *filename);
  void drawCam(int camId);
  void RenderSkeletonToCam(Skeleton *m_pSkeleton,int camId=0);
  void setPoseId(int pId_){poseId = pId_;}  
  void Reset();

private:  
  std::vector<PinholeCamera> allCams;  
  int poseId;
};

#endif

