#include"cameras.h"
#include <FL/gl.h>
#include <FL/glut.H>
#include"skeleton.h"


Cameras::Cameras(){
    poseId = 0;
    allCams.clear();
}
Cameras::~Cameras(){
    allCams.clear();
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



    glColor3f(1,0,0);
    glLineWidth(2.0);
    glBegin(GL_LINES);

    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(0, 0, 0);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);

    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);

    glVertex3f(m_sz * (m_w - 1 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);

    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (m_h - 1 - m_cy) / m_fy, m_sz);
    glVertex3f(m_sz * (0 - m_cx) / m_fx, m_sz * (0 - m_cy) / m_fy, m_sz);

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
    
    

}