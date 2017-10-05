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

#include<Eigen/Core>
#include<Eigen/Dense>





class pinholeCamera{
public:
    pinholeCamera();
    ~pinholeCamera();

private:
    //instrint model:
    float w;
    float h;
    float fx,fy,cx,cy;
    float sz;
    //trajectories:    
};


class Cameras 
{
  //function members
public:    
  Cameras(char *amc_filename, double scale, Skeleton * pSkeleton);

  
  

  ~Cameras();


protected:
  pinholeCamera * allCams;
  int readCameraFile(char* name, double scale);
};

#endif

