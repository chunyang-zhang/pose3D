/*

Revision 1 - Steve Lin (CMU), Jan 14, 2002
Revision 2 - Alla Safonova and Kiran Bhat (CMU), Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao (USC), Feb, 2012

*/
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "skeleton.h"
#include "transform.h"

#ifdef WIN32
#pragma warning(disable : 4996)
#endif

int Skeleton::numBonesInSkel(Bone *bone)
{
  if (bone == NULL)
  {
    return 0;
  }
  else
  {
    return 1 + numBonesInSkel(bone->child) + numBonesInSkel(bone->bros);
  }
}

int Skeleton::movBonesInSkel(Bone *bone)
{
  if (bone == NULL)
  {
    return 0;
  }
  else
  {
    if (bone->dof > 0)
    {
      return 1 + movBonesInSkel(bone->child) + movBonesInSkel(bone->bros);
    }
    else
    {
      return movBonesInSkel(bone->child) + movBonesInSkel(bone->bros);
    }
  }
}

// helper function to convert ASF part name into bone index
int Skeleton::name2idx(std::string name)
{
  int i = 0;
  while (m_pBoneList[i].name != name && i++ < NUM_BONES_IN_ASF_FILE)
    ;
  return m_pBoneList[i].idx;
}

std::string Skeleton::idx2name(int idx)
{
  int i = 0;
  while (m_pBoneList[i].idx != idx && i++ < NUM_BONES_IN_ASF_FILE)
    ;
  return m_pBoneList[i].name;
}

int Skeleton::readASFfile(char *asf_filename, double scale)
{
  //open file
  std::ifstream infile(asf_filename);
  std::string curLine;
  std::stringstream ss;
  std::string keyword;
  if (infile.fail())
    return -1;

  //skip the headers
  while (std::getline(infile, curLine))
  {

    ss.clear();
    ss.str(std::string());
    ss << curLine;
    ss >> keyword;
    if (keyword == ":root")
      break;
  }

  //read the root data
  while (std::getline(infile, curLine))
  {
    std::string keyword;
    ss.str(std::string());
    ss.clear();
    ss << curLine;
    ss >> keyword;
    if (keyword == ":bonedata")
      break;
    if (keyword == "position")
      ss >> this->tx >> this->ty >> this->tz;
    if (keyword == "orientation")
      ss >> this->rx >> this->ry >> this->rz;
  }

  printf("start to get all the orientation\n");
  while (std::getline(infile, curLine))
  {

    std::string keyword;
    ss.str(std::string());
    ss.clear();
    ss << curLine;
    ss >> keyword;
    if (keyword == "begin")
    {
      Bone m_Bone;
      m_Bone.dof = 0;
      for (int dindex = 0; dindex < MAX_DOFS; dindex++)
      {
        m_Bone.dofmask[dindex] = 0;
        m_Bone.dofval[dindex] = 0;
      }
      m_Bone.father = NULL;
      m_Bone.child = NULL;
      m_Bone.bros = NULL;
      m_Bone.name = "";

      while (std::getline(infile, curLine))
      {
        ss.str(std::string());
        ss.clear();
        ss << curLine;
        ss >> keyword;

        if (keyword == "id")
          ss >> m_Bone.idx;
        if (keyword == "name")
          ss >> m_Bone.name;
        if (keyword == "direction")
          ss >> m_Bone.dir[0] >> m_Bone.dir[1] >> m_Bone.dir[2];
        if (keyword == "length")
          ss >> m_Bone.length;
        if (keyword == "axis")
          ss >> m_Bone.axis_x >> m_Bone.axis_y >> m_Bone.axis_z;
        if (keyword == "dof")
        {
          std::string dofstr;
          int tmpdof = 0;
          while (ss >> dofstr)
          {
            m_Bone.dofo[tmpdof] = dofname[dofstr];
            m_Bone.dofmask[dofname[dofstr]] = 1;
            tmpdof++;
            m_Bone.dof++;
          }
        }
        if (keyword == "end")
        {

          m_Bone.dofrx = m_Bone.dofmask[0];
          m_Bone.dofry = m_Bone.dofmask[1];
          m_Bone.dofrz = m_Bone.dofmask[2];
          m_Bone.doftx = m_Bone.dofmask[3];
          m_Bone.dofty = m_Bone.dofmask[4];
          m_Bone.doftz = m_Bone.dofmask[5];
          m_Bone.doftl = m_Bone.dofmask[6];
          m_Bone.length = m_Bone.length * scale;
          m_pBoneList[m_Bone.idx].length = m_Bone.length * scale;
          m_pBoneList[m_Bone.idx] = m_Bone;
          NUM_BONES_IN_ASF_FILE++;
          if (m_Bone.dof > 0)
            MOV_BONES_IN_ASF_FILE++;
          break;
        }
      }
    }
    if (keyword == ":hierarchy")
      break;
  }

  printf("READ %d\n", NUM_BONES_IN_ASF_FILE);

  //Assign parent/child relationship to the bones
  std::getline(infile, curLine);
  ss.str(std::string());
  ss.clear();
  ss << curLine;
  ss >> keyword;
  if (keyword == "begin")
  {
    while (std::getline(infile, curLine))
    {
      ss.str(std::string());
      ss.clear();
      ss << curLine;
      std::string pBoneName, cBoneName;
      int idxP, idxC, idxT;
      ss >> pBoneName >> cBoneName;
      if (pBoneName == "end")
        break;
      idxP = name2idx(pBoneName);
      idxC = name2idx(cBoneName);
      m_pBoneList[idxP].child = &m_pBoneList[idxC];
      m_pBoneList[idxC].father = &m_pBoneList[idxP];
      m_pBoneList[idxC].bros = NULL;
      while (ss >> cBoneName)
      {
        if (cBoneName == "")
          break;
        idxT = name2idx(cBoneName);
        m_pBoneList[idxC].bros = &m_pBoneList[idxT];
        m_pBoneList[idxT].father = &m_pBoneList[idxP];
        m_pBoneList[idxT].bros = NULL;

        idxC = idxT;
      }
    }
  }
  std::cout << "READ END" << std::endl;
  std::cout << "DEBUG: Print tree of Body" << std::endl;

  PrintSkeletonStructure(getRoot());
  //
  std::cout << "Debug: printf Dofo" << std::endl;
  for (int i = 0; i < NUM_BONES_IN_ASF_FILE; i++)
  {
    std::cout << m_pBoneList[i].idx << ":" << m_pBoneList[i].name << ":";
    for (int j = 0; j < m_pBoneList[i].dof; j++)
    {
      switch (m_pBoneList[i].dofo[j])
      {
      case 0:
        std::cout << "rx ";
        break;
      case 1:
        std::cout << "ry ";
        break;
      case 2:
        std::cout << "rz ";
        break;
      case 3:
        std::cout << "tx ";
        break;
      case 4:
        std::cout << "ty ";
        break;
      case 5:
        std::cout << "tz ";
        break;
      case 6:
        std::cout << "tl ";
        break;
      }
    }
    std::cout << std::endl;
  }

  infile.close();
  return 0;
}

/*
This recursive function traverces skeleton hierarchy 
and returns a pointer to the bone with index - bIndex
ptr should be a pointer to the root node 
when this function first called
*/
Bone *Skeleton::getBone(int bIndex)
{
  return &m_pBoneList[bIndex];
}

/*
This function sets father or child for parent bone
If parent bone does not have a child, 
then pChild is set as parent's child
else pChild is set as a father of parents already existing child
Return the pointer to the root bone
*/
Bone *Skeleton::getRoot()
{
  return (m_pRootBone);
}

/***************************************************************************************
Compute relative orientation and translation between the 
parent and child bones. That is, represent the orientation 
matrix and translation vector in the local coordinate of parent body 
*****************************************************************************************/

/*
This function sets rot_parent_current data member.
Rotation from this bone local coordinate system 
to the coordinate system of its parent
*/
void Skeleton::compute_rotation_parent_child(Bone *parent, Bone *child)
{
  if (child != NULL)
  {
    Eigen::Matrix3d parentRot;
    Eigen::Matrix3d childRot;

    parentRot = Eigen::AngleAxisd(-parent->axis_x * M_PI / 180, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-parent->axis_y * M_PI / 180, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(-parent->axis_z * M_PI / 180, Eigen::Vector3d::UnitZ());

    childRot = Eigen::AngleAxisd(child->axis_z * M_PI / 180, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(child->axis_y * M_PI / 180, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(child->axis_x * M_PI / 180, Eigen::Vector3d::UnitX());
    child->rot2parent.setIdentity();
    child->rot2parent.topLeftCorner(3, 3) = parentRot * childRot;
  }
}

// loop through all bones to calculate local coordinate's direction vector and relative orientation
void Skeleton::ComputeBoneRotationToParentCoordSystem()
{
  //compute root to it's parent
  Eigen::Matrix3d rootMat;
  Bone *rootBone = getRoot();
  rootMat = Eigen::AngleAxisd(rootBone->axis_z * M_PI / 180, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(rootBone->axis_y * M_PI / 180, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(rootBone->axis_x * M_PI / 180, Eigen::Vector3d::UnitX());
  rootBone->rot2parent.setIdentity();
  rootBone->rot2parent.topLeftCorner(3, 3) = rootMat;

  //Compute rot_parent_current for all other bones
  int numbones = numBonesInSkel(getRoot());

  for (int i = 0; i < numbones; i++)
  {
    if (m_pBoneList[i].father != NULL)
    {
      compute_rotation_parent_child(m_pBoneList[i].father, &m_pBoneList[i]);
    }
  }
}

/*
Transform the direction vector (dir), 
which is defined in character's global coordinate system in the ASF file, 
to local coordinate
*/
void Skeleton::RotateBoneDirToLocalCoordSystem()
{
  int i;
  for (i = 1; i < NUM_BONES_IN_ASF_FILE; i++)
  {
    Eigen::Vector3d m_BoneDirHomo;
    for (int j = 0; j < 3; j++)
    {
      m_BoneDirHomo(j) = m_pBoneList[i].dir[j];
    }
    Eigen::Matrix3d rotMat;
    rotMat = Eigen::AngleAxisd((-m_pBoneList[i].axis_x) * M_PI / 180, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd((-m_pBoneList[i].axis_y) * M_PI / 180, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd((-m_pBoneList[i].axis_z) * M_PI / 180, Eigen::Vector3d::UnitZ());
    m_BoneDirHomo = rotMat * m_BoneDirHomo;
    for (int j = 0; j < 3; j++)
    {
      m_pBoneList[i].dir[j] = m_BoneDirHomo(j);
    }
  }
}

/******************************************************************************
Interface functions to set the pose of the skeleton 
******************************************************************************/

//Initial posture Root at (0,0,0)
//All bone rotations are set to 0
void Skeleton::setBasePosture()
{
  int i;
  m_RootPos[0] = m_RootPos[1] = m_RootPos[2] = 0.0;

  for (i = 0; i < NUM_BONES_IN_ASF_FILE; i++)
  {
    m_pBoneList[i].rx = m_pBoneList[i].ry = m_pBoneList[i].rz = 0.0;
    m_pBoneList[i].tx = m_pBoneList[i].ty = m_pBoneList[i].tz = 0.0;
  }
  Eigen::Matrix3d iMax;
  iMax.setIdentity();
  calGlobalBonePos(getRoot(), iMax);
}

void Skeleton::enableAllRotationalDOFs()
{
  for (int j = 0; j < NUM_BONES_IN_ASF_FILE; j++)
  {
    if (m_pBoneList[j].dof == 0)
      continue;

    if (!m_pBoneList[j].dofrx)
    {
      m_pBoneList[j].dofrx = 1;
      m_pBoneList[j].rx = 0.0;
      m_pBoneList[j].dof++;
      m_pBoneList[j].dofo[m_pBoneList[j].dof - 1] = 1;
      m_pBoneList[j].dofo[m_pBoneList[j].dof] = 0;
    }

    if (!m_pBoneList[j].dofry)
    {
      m_pBoneList[j].dofry = 1;
      m_pBoneList[j].ry = 0.0;
      m_pBoneList[j].dof++;
      m_pBoneList[j].dofo[m_pBoneList[j].dof - 1] = 2;
      m_pBoneList[j].dofo[m_pBoneList[j].dof] = 0;
    }

    if (!m_pBoneList[j].dofrz)
    {
      m_pBoneList[j].dofrz = 1;
      m_pBoneList[j].rz = 0.0;
      m_pBoneList[j].dof++;
      m_pBoneList[j].dofo[m_pBoneList[j].dof - 1] = 3;
      m_pBoneList[j].dofo[m_pBoneList[j].dof] = 0;
    }
  }
}

// set the skeleton's pose based on the given posture
void Skeleton::setPosture(Posture posture)
{
  m_RootPos[0] = posture.root_pos(0);
  m_RootPos[1] = posture.root_pos(1);
  m_RootPos[2] = posture.root_pos(2);

  for (int j = 0; j < NUM_BONES_IN_ASF_FILE; j++)
  {
    // if the bone has rotational degree of freedom in x direction
    if (m_pBoneList[j].dofrx)
      m_pBoneList[j].rx = posture.bone_rotation[j](0);

    if (m_pBoneList[j].doftx)
      m_pBoneList[j].tx = posture.bone_translation[j](0);

    // if the bone has rotational degree of freedom in y direction
    if (m_pBoneList[j].dofry)
      m_pBoneList[j].ry = posture.bone_rotation[j](1);

    if (m_pBoneList[j].dofty)
      m_pBoneList[j].ty = posture.bone_translation[j](1);

    // if the bone has rotational degree of freedom in z direction
    if (m_pBoneList[j].dofrz)
      m_pBoneList[j].rz = posture.bone_rotation[j](2);

    if (m_pBoneList[j].doftz)
      m_pBoneList[j].tz = posture.bone_translation[j](2);

    if (m_pBoneList[j].doftl)
      m_pBoneList[j].tl = posture.bone_length[j];
  }
  Eigen::Matrix3d iMax;
  iMax.setIdentity();
  calGlobalBonePos(getRoot(), iMax);
}

//Set the aspect ratio of each bone
void Skeleton::set_bone_shape(Bone *bone)
{
  int root = Skeleton::getRootIndex();
  bone[root].aspx = 1;
  bone[root].aspy = 1;
  printf("READ %d\n", numBonesInSkel(bone));
  printf("MOV %d\n", movBonesInSkel(bone));
  int numbones = numBonesInSkel(bone);
  for (int j = 1; j < numbones; j++)
  {
    bone[j].aspx = 0.25;
    bone[j].aspy = 0.25;
  }
}

// Constructor
Skeleton::Skeleton(char *asf_filename, double scale)
{
  dofname.clear();
  dofname.insert(std::pair<std::string, int>("rx", 0));
  dofname.insert(std::pair<std::string, int>("ry", 1));
  dofname.insert(std::pair<std::string, int>("rz", 2));
  dofname.insert(std::pair<std::string, int>("tx", 3));
  dofname.insert(std::pair<std::string, int>("ty", 4));
  dofname.insert(std::pair<std::string, int>("tz", 5));
  dofname.insert(std::pair<std::string, int>("tl", 6));

  NUM_BONES_IN_ASF_FILE = 1;
  MOV_BONES_IN_ASF_FILE = 1;

  m_pBoneList[0].name = "root";
  m_pBoneList[0].idx = 0;
  m_pBoneList[0].dofo[0] = 3;
  m_pBoneList[0].dofo[1] = 4;
  m_pBoneList[0].dofo[2] = 5;
  m_pBoneList[0].dofo[3] = 0;
  m_pBoneList[0].dofo[4] = 1;
  m_pBoneList[0].dofo[5] = 2;
  m_pBoneList[0].dof = 6;
  //Initialization
  m_pBoneList[0].idx = getRootIndex(); // root of hierarchy
  m_pRootBone = &m_pBoneList[0];
  m_pBoneList[0].father = NULL;
  m_pBoneList[0].child = NULL;
  m_pBoneList[0].bros = NULL;
  m_pBoneList[0].dir[0] = 0;
  m_pBoneList[0].dir[1] = 0.;
  m_pBoneList[0].dir[2] = 0.;
  m_pBoneList[0].axis_x = 0;
  m_pBoneList[0].axis_y = 0.;
  m_pBoneList[0].axis_z = 0.;
  m_pBoneList[0].length = 0.05;

  for (int i = 0; i < MAX_DOFS - 1; i++)
  {
    m_pBoneList[0].dofmask[i] = 1;
  }
  m_pBoneList[0].doftx = m_pBoneList[0].dofty = m_pBoneList[0].doftz = 1;
  m_pBoneList[0].dofrx = m_pBoneList[0].dofry = m_pBoneList[0].dofrz = 1;
  m_pBoneList[0].dofmask[MAX_DOFS - 1] = 0;
  m_RootPos[0] = m_RootPos[1] = m_RootPos[2] = 0;
  //	m_NumDOFs=6;
  tx = ty = tz = rx = ry = rz = 0.0;
  // build hierarchy and read in each bone's DOF information
  int code = readASFfile(asf_filename, scale);
  if (code != 0)
    throw 1;

  //transform the direction vector for each bone from the world coordinate system
  //to it's local coordinate system

  std::cout << "DEBUG: Start Rotate Bone Dir to Local" << std::endl;
  RotateBoneDirToLocalCoordSystem();
  std::cout << "DEBUG: End Rotate Bone Dir to Local" << std::endl;

  //Calculate rotation from each bone local coordinate system to the coordinate system of its parent
  //store it in rot_parent_current variable for each bone

  ComputeBoneRotationToParentCoordSystem();

  std::cout << "DEBUG: End Rotate Bone Dir to Parent" << std::endl;

  //

  //Set the aspect ratio of each bone
  set_bone_shape(m_pRootBone);
}

Skeleton::~Skeleton()
{
}

void Skeleton::GetRootPosGlobal(double rootPosGlobal[3])
{
  rootPosGlobal[0] = m_RootPos[0];
  rootPosGlobal[1] = m_RootPos[1];
  rootPosGlobal[2] = m_RootPos[2];
}

void Skeleton::GetTranslation(double translation[3])
{
  translation[0] = tx;
  translation[1] = ty;
  translation[2] = tz;
}

void Skeleton::GetRotationAngle(double rotationAngle[3])
{
  rotationAngle[0] = rx;
  rotationAngle[1] = ry;
  rotationAngle[2] = rz;
}
void Skeleton::PrintSkeletonStructure(Bone *bone, int indent)
{
  if (bone != NULL)
  {
    if (indent)
      std::cout << std::setw(indent) << " ";
    std::cout << bone->name << "\n";
    PrintSkeletonStructure(bone->child, indent + 4);
    PrintSkeletonStructure(bone->bros, indent);
  }
}

void Skeleton::calGlobalBonePos(Bone *m_bone, Eigen::Matrix3d f2w)
{
  if (m_bone == NULL)
    return;

  calGlobalBonePos(m_bone->bros, f2w);

  Eigen::Vector3d jointPoseInLocal = Eigen::Map<Eigen::Vector3d>(m_bone->dir)* m_bone->length;

  double crz = 0.0;
  double cry = 0.0;
  double crx = 0.0;

  if (m_bone->dofrz)
    crz = m_bone->rz;
  if (m_bone->dofry)
    cry = m_bone->ry;
  if (m_bone->dofrx)
    crx = m_bone->rx;

  Eigen::Matrix3d c2a;
  c2a = Eigen::AngleAxisd(crz * M_PI / 180, Eigen::Vector3d::UnitZ()) 
    * Eigen::AngleAxisd(cry * M_PI / 180, Eigen::Vector3d::UnitY())
    * Eigen::AngleAxisd(crx * M_PI / 180, Eigen::Vector3d::UnitX());
  
  Eigen::Matrix3d a2f = m_bone->rot2parent.topLeftCorner(3, 3);
  Eigen::Matrix3d c2w = c2a * a2f * f2w;



  
  


  Eigen::Vector3d fatherPos; 
  if (m_bone->idx == 0){
    fatherPos = Eigen::Map<Eigen::Vector3d>(m_RootPos); 
  }
  else{
    fatherPos = (m_bone->father)->jointPos;
  } 
  
  

  m_bone->jointPos = f2w*a2f*c2a*jointPoseInLocal + fatherPos;                
  
  
  calGlobalBonePos(m_bone->child, f2w*a2f*c2a);    
}