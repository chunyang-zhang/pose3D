/*

Revision 1 - Steve Lin (CMU), Jan 14, 2002
Revision 2 - Alla Safonova and Kiran Bhat (CMU), Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao (USC), Feb, 2012

*/
#include <stdio.h>
#include <string.h>
#include <fstream>
#include<sstream>
#include <math.h>
#include <stdlib.h>


#include "skeleton.h"
#include "motion.h"

Motion::Motion(char *amc_filename, double scale, Skeleton * pSkeleton_)
{
  pSkeleton = pSkeleton_;
  m_NumFrames = 0;
  m_pPostures.clear();
  int code = readAMCfile(amc_filename, scale);	
  if (code < 0)
    throw 1;
}

Motion::~Motion()
{
  m_pPostures.clear();
}



//Set posture at spesified frame
void Motion::SetPosture(int frameIndex, Posture InPosture)
{
  m_pPostures[frameIndex] = InPosture; 	
}

void Motion::SetBoneRotation(int frameIndex, int boneIndex, Eigen::Vector3d vRot)
{
  m_pPostures[frameIndex].bone_rotation[boneIndex] = vRot;
}

void Motion::SetRootPos(int frameIndex, Eigen::Vector3d vPos)
{
  m_pPostures[frameIndex].root_pos = vPos;
}

Posture * Motion::GetPosture(int frameIndex)
{
  if (frameIndex < 0 || frameIndex >= m_NumFrames)
  {
    printf("Error in Motion::GetPosture: frame index %d is illegal.\n", frameIndex);
    printf("m_NumFrames = %d\n", m_NumFrames);
    exit(0);
  }
  return &(m_pPostures[frameIndex]);
}

int Motion::readAMCfile(char* name, double scale)
{
  Bone *hroot, *bone;
  bone = hroot = pSkeleton->getRoot();


  std::string curLine;
  std::string bname;
  std::stringstream ss;
  std::ifstream infile(name, std::ios::in );  
  if(infile.fail() ) return -1;
          
  //Compute number of frames. 
  //Subtract 3 to  ignore the header
  //There are (NUM_BONES_IN_ASF_FILE - 2) moving bones and 2 dummy bones (lhipjoint and rhipjoint)
  int numbones = pSkeleton->numBonesInSkel(bone);
  int movbones = pSkeleton->movBonesInSkel(bone);

  //skip the header:  
  getline(infile,curLine);
  //reserved to be used later
  getline(infile,curLine);
  getline(infile,curLine);

  int nFrames=0;

  while(getline(infile,curLine)){        
    Posture m_post;
    ss.clear();        
    ss.str(std::string());
    ss<<curLine;
    ss>>m_post.idx;
    
    for(int j=0;j<movbones;j++){
      getline(infile,curLine);
      ss.str(std::string());
      ss.clear();
      ss<<curLine;
      ss>>bname;
      int bone_idx = pSkeleton->name2idx(bname);
      Bone * m_Bone = pSkeleton->getBone(bone_idx);
      m_post.bone_rotation[bone_idx] = Eigen::Vector3d(0,0,0);            
      for(int dofid = 0 ; dofid<m_Bone->dof;dofid++){  
        double tmp;
        ss>>tmp;
        switch(m_Bone->dofo[dofid]){          
        case 0:
          m_post.bone_rotation[bone_idx](0)=tmp;          
          break;
        case 1:
          m_post.bone_rotation[bone_idx](1)=tmp;
          break;
        case 2: 
          m_post.bone_rotation[bone_idx](2)=tmp;
          break;
        case 3:
          m_post.bone_translation[bone_idx](0)=tmp*scale;
          break;
        case 4:
          m_post.bone_translation[bone_idx](1)=tmp*scale;
          break;
        case 5:
          m_post.bone_translation[bone_idx](2)=tmp*scale;
          break;
        case 6:
          m_post.bone_length[bone_idx]=tmp;// * scale;
          break;
        }
      }
    }
    m_post.root_pos(0)= m_post.bone_translation[0](0);
    m_post.root_pos(1)= m_post.bone_translation[0](1);
    m_post.root_pos(2)= m_post.bone_translation[0](2);

    m_pPostures.push_back(m_post);
    nFrames++;    
  }
  infile.close();
  printf("%d samples in '%s' are read.\n", nFrames, name);
  m_NumFrames = nFrames;
  return nFrames;
}