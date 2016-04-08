#ifndef VOXEL_HPP
#define VOXEL_HPP
#include "Mesh.hpp"
#include <vector>
#include <iostream>
#define VOXEL_DIM 3

struct GridIdx{
  int i[3];
  GridIdx(int a=0,int b=0,int c=0){
    i[0]=a;
    i[1]=b;
    i[2]=c;
  }
  GridIdx&operator=(const GridIdx & ii){
    for(int jj =0;jj<3;jj++){
      i[0] = ii.i[0];
    }
    return *this;
  }

  int & operator[](int idx){return i[idx];}
  int operator[](int idx)const{return i[idx];}

  bool operator<(const GridIdx & grid)const{
    for(int ii=0;ii<VOXEL_DIM;ii++){
      if(i[ii]==grid.i[ii]){
        continue;
      }
      if(i[ii]<grid.i[ii]){
        return true;
      }
      return false;
    }
    return false;
  }
};

struct Grid
{
  int * fill;
  Grid():fill(0),
    nx(0),ny(0),nz(0){}
  ~Grid(){
    if(fill != 0){
      delete [] fill;
    }
  }
  void allocate(int _nx,int _ny,int _nz){
    if(fill != 0){
      delete [] fill;
    }
    nx = _nx;
    ny = _ny;
    nz = _nz;
    int size = nx*ny*nz;
    fill = new int[size];
    for(int ii = 0;ii<size;ii++){
      fill[ii] = 0;
    }
  }

  int & operator()(int ii, int jj, int kk){
    if(ii>=nx || ii<0 ||
       jj>=ny || jj<0 ||
       kk>=nz || kk<0){
      std::cout<<"Error: Index out of bound "<<ii<<" "<<jj<<" "<<kk<<"\n";
    }
    return fill[ii*ny*nz + jj*nz + kk];
  }

  int read(int ii, int jj, int kk)const{
    if(ii>=nx || ii<0 ||
       jj>=ny || jj<0 ||
       kk>=nz || kk<0){
      return 0;
    }
    return fill[ii*ny*nz + jj*nz + kk];
  }

  int nx,ny,nz;
  int size(int dim){
    switch (dim){
    case 0:
      return nx;
    case 1:
      return ny;
    case 2:
      return nz;
    }
    return 0;
  }
};

///@brief builds a voxel grid for a mesh where 0 stands for empty cell
///and 1 stands for a filledd cell.
class Voxel{
public:

  typedef Vector3f Vec3f;

  Voxel();
  ~Voxel();

  Grid grid;
  float gridlen;
  Vec3f orig;
  int res;

  int nx, ny, nz;

  ///@brief additional layer of voxels around the bounding box of the object
  int voxelPadding;

  void buildGrid();

  ///@brief allocates arrays and sets origin
  void initGrid(Vector3f mn, Vector3f mx);

  ///@brief insert a triangle
  ///@param tidx triangle index
  void insert(int tidx, Grid & inside);

  ///@brief flood fill the exterior to figure out the interior
  void floodfill(Grid & inside);
  void saveObj(const char * outfile);
  void computeConnectivity();

  ///@brief save a 3D voxel mesh for simulation.
  void saveConnectivity(const char * outfile,
                        bool saveParts=false, bool reflectxz=false);

  ///@brief slice in Z direction out of the screen.
  void save2DSlice(const char * outfile, int slice=0,
                   bool saveParts = false);

  Mesh * m;

private:
  static Mesh cube;

  int clampIdx(int idx, int dim);
  void bbox(int ii, int * tmn, int * tmx);

  void rasterize(int tidx, Grid & grid, Grid & inside);

  ///@brief computes grid index given coordinate
  void vec2grid(const Vec3f &v, GridIdx & grid);
  bool trigCubeIntersect(int tidx, const GridIdx & cube);
  int  trigCubeSign     (int tidx, const GridIdx & cube);
  //maximum grid size
  Vec3i mg;
  
  ///@brief vertices of the grid that are part of the hex element mesh.
  std::vector<Vector3f> ev;

  ///@brief index of a vertex at (i,j,k) to the vertex list v
  ///-1 if the vertex is not part of the mesh.
  std::vector<int> vidx;
  std::vector<std::vector< int> > e;

  ///@brief For a voxel grid with multiple parts
  std::vector<int> partIdx;

  ///@brief index into e for a grid point EI(i,j,k).
  std::vector<int> gridToE;

  int EI(int ii, int jj, int kk){
    return grid.ny * grid.nz * ii + grid.nz * jj + kk;
  }

  int VI(int ii, int jj, int kk){
    return (grid.ny+1) * (grid.nz + 1 ) * ii + 
      (grid.nz + 1) * jj + kk;
  }

};
#endif
