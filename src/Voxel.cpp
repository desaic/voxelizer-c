#include "Voxel.hpp"
#include "trigAABBIntersect.hpp"
#include <queue>
#include <set>
#include <cstdlib>
#include <fstream>
Mesh Voxel::cube;
extern Mesh UNIT_CUBE;
const int n_nbr= 6;
int nbr[n_nbr][3] = {{-1, 0, 0},{1,0,0},
                 { 0,-1, 0},{0,1,0},
                 { 0, 0,-1},{0,0,1}};

int QUAD_V[4][2]={
  {0,0},
  {0,1},
  {1,0},
  {1,1}
};

///@brief vertex ordering for one hex element.
int eVidx[8][3]={
{0,0,0},
{0,0,1},
{0,1,0},
{0,1,1},
{1,0,0},
{1,0,1},
{1,1,0},
{1,1,1}
};

void BBox(const std::vector<Vector3f >& v,
    Vector3f & mn, Vector3f & mx);

bool inbound(int x, int lb, int ub)
{
  return x>=lb && x<=ub;
}

void
Voxel::saveObj(const char * outfile)
{
  Mesh box;
  Mesh mout;
  Vector3f dx(1,1,1);
  int nx = grid.nx;
  int ny = grid.ny;
  int nz = grid.nz;
  for (int ii = 0; ii < nx; ii++) {
  for (int jj = 0; jj < ny; jj++) {
  for (int kk = 0; kk < nz; kk++) {
    if(grid(ii,jj,kk)==0){
      continue;
    }
    Vector3f coord(0.5f + ii, 0.5f + jj, 0.5f+kk);
    Vector3f box0 = gridlen * (coord - 0.5f * dx) + orig;
    Vector3f box1 = gridlen * (coord + 0.5f * dx) + orig;
    makeCube(box, box0, box1);
    mout.append(box);
  }
  }
  }
  mout.save_obj(outfile);
}

void Voxel::computeConnectivity()
{
  int nvx = grid.nx+1;
  int nvy = grid.ny+1;
  int nvz = grid.nz+1;
  vidx.resize(nvx*nvy*nvz, -1);
  gridToE.resize(grid.nx*grid.ny*grid.nz,-1);
  for (int ii = 0; ii < grid.nx; ii++) {
  for (int jj = 0; jj < grid.ny; jj++) {
  for (int kk = 0; kk < grid.nz; kk++) {
    if(grid(ii,jj,kk)==0){
      continue;
    }
    std::vector<int> nodes(8,0);
    for(int ll = 0; ll<8; ll++){
      int vi = ii + eVidx[ll][0];
      int vj = jj + eVidx[ll][1];
      int vk = kk + eVidx[ll][2];
      //index of vertex in list of grid vertices
      int gridv = VI(vi,vj,vk);
      //index of vertex in the hex mesh
      int meshv = vidx[gridv];
      if(meshv<0){
        vidx[gridv] = ev.size();
        meshv = vidx[gridv];
        Vector3f coord = orig + gridlen * Vector3f(vi,vj,vk);
        ev.push_back(coord);
      }
      
      nodes[ll] = meshv;
    }
    int eIdx = EI(ii,jj,kk);
    int p = grid(ii,jj,kk);
//    std::cout<<e.size()<<"\n";
    gridToE[eIdx] = e.size();
    e.push_back(nodes);
    partIdx.push_back(p);
  }
  }
  }
}

void Voxel::saveConnectivity(const char * outfile, bool saveParts,
                             bool reflectxz)
{
  std::ofstream out(outfile);
  if(!out.good()){
    std::cout<<"Cannot open "<<outfile<<"\nin saveConnectivity\n";
    return;
  }
  out<<"#verts "<<ev.size()<<"\n";
  out<<"#elts " << e.size()<<"\n";
  float factor = 1;
  if(reflectxz){
    factor = -1;
  }
  for(unsigned int ii =0; ii<ev.size(); ii++){
    out<< factor*ev[ii][0] <<" "<< ev[ii][1] << " " << factor*ev[ii][2] <<"\n";
  }
  for(unsigned int ii =0 ; ii< e.size(); ii++){
    out<<"8";
    for(int jj =0 ; jj<8; jj++){
      out<<" "<<e[ii][jj];
    }
    out<<"\n";
  }

  if(saveParts){
    out<<"parts\n";
    out<<partIdx.size()<<"\n";
    for(unsigned int ii = 0; ii<partIdx.size(); ii++){
      out<<partIdx[ii]<<"\n";
    }
    out<<"density\n";
    out<<e.size()<<"\n";
    //trivial for 3d
    for(unsigned int ii =0; ii<e.size(); ii++){
      out<<"1\n";
    }
  }
  out.close();
}

void Voxel::save2DSlice(const char * outfile, int slice,
                        bool saveParts)
{
  std::ofstream out(outfile);
  if(!out.good()){
    std::cout<<"Cannot open "<<outfile<<"\nin saveConnectivity\n";
    return;
  }

  const int N_QUADV=4;
  //index of a vertex in the list of vertices used by the slice.
  std::vector<std::vector<int> > quadVidx;
  std::vector<std::vector<int> > quads;
  std::vector<Vector3f> quadv;
  std::vector<int> quadPart;
  std::vector<float> quadDensity;
  int vCnt = 0;

  quadVidx.resize(grid.nx + 1);
  for (int ii = 0; ii<=grid.nx; ii++) {
    quadVidx[ii].resize(grid.ny + 1, -1);
  }

  for (int ii = 0; ii<grid.nx; ii++) {
  for (int jj = 0; jj < grid.ny; jj++) {
    bool filled = false;
    float density = 0;
    int eidx = -1;
    for(int kk = 0; kk< grid.nz; kk++){
      if(grid(ii,jj,kk)!=0){
        filled = true;
        eidx = gridToE[EI(ii,jj,kk)];
        density+=1.0;
      }
    }
    if(!filled){
      continue;
    }

    std::vector<int> q(N_QUADV,0);
    for(int vv = 0; vv<N_QUADV; vv++){
      int vi = ii + QUAD_V[vv][0];
      int vj = jj + QUAD_V[vv][1];
//      std::cout<<vidx<<"\n";
      int vidx = quadVidx[vi][vj];
      if(vidx<0){
        vidx = vCnt;
        quadVidx[vi][vj] = vidx;
        Vector3f v(vi*gridlen, vj*gridlen, 0);
        quadv.push_back(v);
        vCnt++;
      }
      q[vv] = vidx;
    }
    quads.push_back(q);
    quadPart.push_back(partIdx[eidx]);
    density = density * gridlen;
    quadDensity.push_back(density);
  }
  }

  out<<"#verts "<<vCnt<<"\n";
  out<<"#quads "<<quads.size()<<"\n";
  for(unsigned int ii =0; ii<quadv.size(); ii++){
    out<< quadv[ii][0] <<" "<<quadv[ii][1]<<" "<<0<<"\n";
  }

  for(unsigned int ii =0 ; ii<quads.size(); ii++){
    out<<quads[ii].size();
    for(unsigned int jj =0 ; jj<quads[ii].size() ; jj++){
      out<<" "<<quads[ii][jj];
    }
    out<<"\n";
  }

  if(saveParts){
    out<<"parts\n";
    out<<quadPart.size()<<"\n";
    for(unsigned int ii = 0; ii<quadPart.size(); ii++){
      out<<quadPart[ii]<<"\n";
    }

    out<<"density\n";
    out<<quadDensity.size()<<"\n";
    for(unsigned int ii = 0; ii<quadDensity.size(); ii++){
      out<<quadDensity[ii]<<"\n";
    }
  }

}

void Voxel::floodfill(Grid & inside)
{
  //find an empty cell
  int ii=0,jj=0,kk=0;
  for(ii= 0;ii<grid.nx; ii++){
    for(jj= 0; jj<grid.ny; jj++){
      for(kk= 0; kk<grid.nz; kk++){
        if(grid(ii,jj,kk)==0){
          goto found_empty;
        }
      }
    }
  }
  found_empty:
  std::cout << "floodfill found empty \n";
  std::queue<GridIdx> que;
  que.push(GridIdx(ii,jj,kk));
  Grid tmpGrid;
  tmpGrid.allocate(grid.nx,grid.ny,grid.nz);
  tmpGrid(ii,jj,kk) = 1;

  while(!que.empty()){
    GridIdx idx = que.front();
    que.pop();
    for(int ii = 0;ii < n_nbr;ii++){
      GridIdx ni = idx;
      for(int jj = 0;jj<3;jj++){
        ni[jj] += nbr[ii][jj];
      }
      if(!inbound(ni[0],0,grid.nx-1)
       ||!inbound(ni[1],0,grid.ny-1)
       ||!inbound(ni[2],0,grid.nz-1)){
        continue;
      }
      if( (tmpGrid(ni[0],ni[1],ni[2])==0) &&
         (grid(ni[0],ni[1],ni[2])==0)){
        que.push(ni);
        tmpGrid(ni[0],ni[1],ni[2]) = 1;
      }
    }
  }
  for(ii= 0;ii<grid.nx; ii++){
    for( jj= 0; jj<grid.ny; jj++){
      for(kk= 0; kk<grid.nz; kk++){
        if(tmpGrid(ii,jj,kk)==0){
          //if(grid(ii,jj,kk)==0){
          //  inside(ii,jj,kk)=1;
          //}
          grid(ii,jj,kk)=1;
        }
      }
    }
  }
}

void Voxel::insert(int tidx, Grid & inside)
{
  rasterize(tidx, grid, inside);
}

Voxel::~Voxel()
{}

int
maxIdxVec3f(const Voxel::Vec3f & vec){
  float val = vec[0];
  int idx=0;
  for(int ii = 1;ii<3;ii++){
    if(val<vec[ii]){
      idx = ii;
      val = vec[ii];
    }
  }
  return idx;
}

Voxel::Voxel():
  gridlen(-1),
  res(0),nx(0), ny(0), nz(0),
  voxelPadding(2), m(0)
{
  mg[0] = 10000;
  mg[1] = 10000;
  mg[2] = 10000;
}

void
Voxel::initGrid(Vector3f mn, Vector3f mx)
{
  Voxel::cube=UNIT_CUBE;
  Vec3f size = mx-mn;
  int mi = maxIdxVec3f(size);

  //not initialized
  if(gridlen<0){
    gridlen = size[mi]/res;
  }

  nx = size[0]/gridlen + 2*voxelPadding;
  ny = size[1]/gridlen + 2*voxelPadding;
  nz = size[2]/gridlen + 2*voxelPadding;
  nx = std::min(nx,mg[0]);
  ny = std::min(ny,mg[1]);
  nz = std::min(nz,mg[2]);

  grid.allocate(nx,ny,nz);
  orig = mn - (voxelPadding - 1e-4) * Vector3f(gridlen,gridlen,gridlen);

}

void
Voxel::buildGrid()
{
  Grid inside;
  inside.allocate(nx,ny,nz);
  for(size_t ii = 0;ii<m->t.size();ii++){
    //std::cout << ii << "\n";
    insert(ii, inside);
  }
  floodfill(inside);
//  for(int ii =0 ; ii<nx; ii++){
//    for(int jj =0 ; jj<ny; jj++){
//      for(int kk =0 ; kk<nz; kk++){
//        if( (!inside(ii,jj,kk)) && grid(ii,jj,kk)){
//          grid(ii,jj,kk)=false;
//        }
//      }
//    }
//  }
}

///@return 2 if cube center is on the positive half space (outside) of triangle.
///@return 1 of tje cube center is on the negative side, but the projection is outside
///of the triangle.
int Voxel::trigCubeSign(int tidx, const GridIdx & gi)
{
//   std::cout<<gi[0]<<" "<<gi[1]<<" "<<gi[2]<<"\n";
  Vector3f boxcenter((float)((0.5+gi[0])*gridlen+orig[0]),
                      (float)((0.5+gi[1])*gridlen+orig[1]),
                      (float)((0.5+gi[2])*gridlen+orig[2]));
  Vector3f tv[3];
  for(int ii=0;ii<3;ii++){
    tv[ii]=m->v[m->t[tidx][ii]];
  }
  Vector3f e1 = tv[1] - tv[0];
  Vector3f e2 = tv[2] - tv[0];
  Vector3f n = Vector3f::cross(e1,e2).normalized();
  Vector3f d = boxcenter - tv[0];
  float dotp = Vector3f::dot(n,d);
  if(dotp > 0){
    return 2;
  }
  //too far away
  if(dotp < -0.9 * gridlen){
    return 3;
  }
  n.normalize();
  d = boxcenter - (Vector3f::dot(n,d))*n;
  Vector3f n0 = Vector3f::cross(tv[1] - d, tv[2] - d);
//  n0.normalize();
  float thresh = 1;
//for jumpers
//  float thresh = 0.00005;

  if(Vector3f::dot(n0,n)<-thresh){
    return 1;
  }
  
  Vector3f n1 = Vector3f::cross(tv[2] - d, tv[0] - d);
//  n1.normalize();
  if(Vector3f::dot(n1,n)<-thresh){
    return 1;
  }
  
  Vector3f n2 = Vector3f::cross(tv[0] - d, tv[1] - d);
//  n2.normalize();
  if(Vector3f::dot(n2,n)<-thresh){
    return 1;
  }
  
  return -1;
}

bool Voxel::trigCubeIntersect(int tidx, const GridIdx & cube)
{
  float boxcenter[3]={(float)((0.5+cube[0])*gridlen+orig[0]),
                      (float)((0.5+cube[1])*gridlen+orig[1]),
                      (float)((0.5+cube[2])*gridlen+orig[2])};
  float boxhalfsize[3]={(float)gridlen*.5,
                      (float)gridlen*.5,
                      (float)gridlen*.5};
  float triverts[3][3];
  for(int ii=0;ii<3;ii++){
    for(int jj=0;jj<3;jj++){
      triverts[ii][jj]=m->v[m->t[tidx][ii]][jj];
    }
  }
  return triBoxOverlap(boxcenter,boxhalfsize, triverts);
}

int Voxel::clampIdx(int idx, int dim)
{
  idx = std::max(0,idx);
  idx = std::min(idx, grid.size(dim) - 1);
  return idx;
}

void Voxel::vec2grid(const Vec3f & v,  GridIdx & grid)
{
  for(int ii=0; ii<VOXEL_DIM; ii++) {
    grid[ii]= (int)((v[ii]-orig[ii])/gridlen);
    grid[ii] = clampIdx(grid[ii],ii);
  }
}

void Voxel::rasterize(int tidx, Grid & grid, Grid & inside)
{
  //bounding box of a triangle
  int tmn[3], tmx[3];
  bbox(tidx,tmn,tmx);
  for(int ix=tmn[0];ix<=(tmx[0]);ix++){
    for(int iy=tmn[1];iy<=(tmx[1]);iy++){
      for(int iz=tmn[2];iz<=(tmx[2]);iz++){
        GridIdx gi(ix,iy,iz);
        if(trigCubeIntersect(tidx,gi)){
          //check if the center of the voxel is in the positive half plane of the triangle
          //if(trigCubeSign(tidx,gi)<0){
            inside(ix,iy,iz)=true;
            grid(ix,iy,iz)=true;
          //}
        }
      }
    }
  }
}

void Voxel::bbox(int ii, int * tmn, int * tmx)
{
  GridIdx vidx;
  vec2grid(m->v[m->t[ii][0]],vidx);
  for(int jj=0; jj<3; jj++) {
    tmn[jj]=vidx[jj]-1;
    tmx[jj]=vidx[jj];
  }
  for(int jj=1; jj<3; jj++) {
    vec2grid(m->v[m->t[ii][jj]],vidx);
    for(int kk=0; kk<VOXEL_DIM; kk++) {
      if(vidx[kk]-1<tmn[kk]) {
        tmn[kk]=vidx[kk]-1;
      }
      if(vidx[kk]>tmx[kk]) {
        tmx[kk]=vidx[kk];
      }
    }
  }
  for(int ii = 0;ii<VOXEL_DIM;ii++){
    tmn[ii] = std::max(0,tmn[ii]-1);
    tmx[ii] = std::min(grid.size(ii)-2,tmx[ii]+2);
  }
}
