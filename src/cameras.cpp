#include"cameras.h"
#include <FL/gl.h>
#include <FL/glut.H>
#include"skeleton.h"
#include<opencv2/opencv.hpp>
#include<stdlib.h>



Cameras::Cameras(){
    poseId = 0;
    allCams.clear();
    m_BoneColor.clear();
}
Cameras::~Cameras(){
    allCams.clear();
    m_BoneColor.clear();

}
void Cameras::Reset(){
    for(int i=0;i<allCams.size();i++){
        allCams[i].pos.clear();        
    }    
    allCams.clear();
}
int Cameras::loadNewCameraFromFile(char *filename){

    PinholeCamera m_Cam;
    m_Cam.pos.clear();
    m_Cam.poseId = 0;
    std::ifstream infile(filename);
    std::string curLine;
    std::string keyword;
    std::stringstream ss;
    

    //skip the header
    while(std::getline(infile,curLine)){
        ss.clear();
        ss.str(std::string());
        ss<<curLine;
        ss>>keyword;
        if(keyword==":camera") break;        
    }
    // read camera instinct
    double tfx,tfy,tcx,tcy;
    double tw,th,tsz;
    m_Cam.K.setZero();
    std::getline(infile,curLine);
    ss.clear();
    ss.str(std::string());
    ss<<curLine;
    ss>>tw>>th>>tfx>>tfy>>tcx>>tcy>>tsz;
    m_Cam.K(0,0)=tfx;
    m_Cam.K(1,1)=tfy;
    m_Cam.K(0,2)=tcx;
    m_Cam.K(1,2)=tcy;
    m_Cam.K(2,2)=1;
    m_Cam.w = tw;
    m_Cam.h = th;
    m_Cam.sz = tsz;
    // read camera extinct
    while(std::getline(infile,curLine)){
        ss.clear();
        ss.str(std::string());
        ss<<curLine;
        double m_PosArr[16];
        for(int j=0;j<16;j++){
            ss>>m_PosArr[j];
        }
        Eigen::Matrix4d m_Pos = Eigen::Map<Eigen::Matrix4d>(m_PosArr);
        m_Cam.pos.push_back(m_Pos);        
    }
    allCams.push_back(m_Cam);        
    return 0;


}
void Cameras::drawCam(int camId){
    if(!allCams.size())
     return;    
    glPushMatrix();
    GLint lightingStatus;
    glGetIntegerv(GL_LIGHTING, &lightingStatus);
    glDisable(GL_LIGHTING);
  
    
    PinholeCamera mCam_cur = allCams[camId];    
    glMultMatrixd((double *)mCam_cur.pos[poseId].data());

    double m_w,m_h;
    double m_sz;
    double m_fx,m_fy,m_cx,m_cy;
    m_w = mCam_cur.w;    
    m_h = mCam_cur.h;
    m_sz = mCam_cur.sz;
    m_fx = mCam_cur.K(0,0);
    m_fy = mCam_cur.K(1,1);
    m_cx = mCam_cur.K(0,2);
    m_cy = mCam_cur.K(1,2);
        
    glLineWidth(2.0);
    glBegin(GL_LINES);

    glColor3f(0.2,0.2,0.2);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);


    glColor3f(0.2,0.2,0.2);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);

    glColor3f(0.2,0.2,0.2);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);

    glColor3f(0.2,0.2,1.0);    
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz+m_sz);

    glColor3f(0.2,1.0,0.2);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);

    glColor3f(1.2,0.2,0.2);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glEnd();

    if (lightingStatus)
    glEnable(GL_LIGHTING);
    glPopMatrix();
}
void Cameras::RenderSkeletonToCam(Skeleton *m_p,int camId){            
    if(!allCams.size())
        return;    
    //perspective 
    PinholeCamera  m_pCam = allCams[camId];

    //1 camera Pose
    cv::Mat img2show((int)m_pCam.h,(int)m_pCam.w,CV_8UC(3),cv::Scalar::all(55));
    cv::Mat img2showFlip;
    int numJoints = m_p->numBonesInSkel(m_p->getRoot());
    if(!m_BoneColor.size()){
        for(int i=0;i<numJoints;i++){
            cv::Scalar m_color = cv::Scalar(rand()%255,rand()%255,rand()%255);
            m_BoneColor.push_back(m_color);            
        }                
    }
    double * imageX = new double[numJoints];
    double * imageY = new double[numJoints];

    Eigen::Matrix4d cam_pos = m_pCam.pos[poseId];    

    //2 Joint  Pose
    for(int i=0;i<numJoints;i++){
        Eigen::Vector4d jPos;
        jPos(3)=1;
        jPos.head(3) = m_p->getBone(i)->jointPos;
        Eigen::Vector3d jPosInFrame = (cam_pos.inverse()*jPos).head(3);
        jPosInFrame = m_pCam.K*jPosInFrame;
        imageX[i] = jPosInFrame(0)/jPosInFrame(2);
        imageY[i] = jPosInFrame(1)/jPosInFrame(2);        
    }

    for(int i=1;i<numJoints;i++){
        Bone *m_bone = m_p->getBone(i);
        int iref = m_bone->father->idx;
        cv::line(img2show,cv::Point(imageX[i],imageY[i]),cv::Point(imageX[iref],imageY[iref]),m_BoneColor[i],5);
    }    
    cv::flip(img2show,img2showFlip,-1);
    cv::imshow("images",img2showFlip);
    cv::waitKey(1);


    delete [] imageX;
    delete [] imageY;
}