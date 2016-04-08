#include "Voxel.hpp"

#include <stdlib.h>
#include <iostream>
#include <fstream>

int main(int argc ,char * argv[])
{
  int voxRes = 160;
  int nMesh = 1;
  int slice = 6;
  std::vector<std::string> meshfiles;

  if(argc<3){
    std::cout<<argv[0]<<"#meshes meshfiles(obj format) [-r voxelres:int | -l gridlen:float]\n";
    return 0;
  }

  nMesh = std::atoi(argv[1]);
  if(argc<nMesh+2){
    std::cout<<"not enough mesh files\n";
    return 0;
  }

  int argIdx = 2;
  for(int ii = 0; ii<nMesh; ii++){
    meshfiles.push_back(std::string(argv[argIdx]));
    argIdx ++ ;
  }

  bool specifyLen = false;
  float gridlen = 0;
  if(argc>argIdx){
    //specify length
    if(argv[argIdx][1]=='l'){
      argIdx ++;
      specifyLen = true;
      gridlen = atof(argv[argIdx]);
    }else if(argv[argIdx][1] == 'r'){
      argIdx ++;
      //specify number of voxels in longest axis
      voxRes = atoi(argv[argIdx]);
    }
    argIdx ++;
  }

  int dim = 3;
  if(argc>argIdx){
    dim = atoi(argv[argIdx]);
    argIdx ++;
    if(dim==2 && argc>argIdx){
      slice = atoi(argv[argIdx]);
      argIdx ++;
    }
  }
  std::cout<<"Voxel grid resolution: "<<voxRes<<"\n";
  std::vector<Mesh> m(nMesh);

  Vector3f mx, mn;
  for(unsigned int ii = 0; ii<m.size(); ii++){
    std::cout << "Load mesh " << meshfiles[ii].c_str() << "\n";
    m[ii].load_mesh(meshfiles[ii].c_str(), false);
    if(ii==0){
      BBox(m[ii].v, mn, mx);
    }else{
      BBox(m[ii].v, mn, mx, true);
    }
  }

  Voxel voxelunion;
  if(specifyLen){
    voxelunion.gridlen = gridlen;
  }
  voxelunion.res = voxRes;
  voxelunion.initGrid(mn, mx);
  Vector3f bsize = mx - mn;
  mx = mx + 1e-5 * bsize;
  mn = mn - 1e-5 * bsize;
  std::cout << bsize[0] << " size\n";
  int vx = voxelunion.grid.nx;
  int vy = voxelunion.grid.ny;
  int vz = voxelunion.grid.nz;

  for(int mm =0 ; mm<nMesh; mm++){
    Voxel voxel;
    voxel.res = voxRes;
    if(specifyLen){
      voxel.gridlen = gridlen;
    }
    voxel.initGrid(mn, mx);
    voxel.m = &m[mm];
    voxel.buildGrid();
    for(int ii = 0;ii<vx;ii++){
      for(int jj = 0;jj<vy;jj++){
        for(int kk = 0;kk<vz;kk++){
          if(voxel.grid.read(ii,jj,kk)>0 && voxelunion.grid(ii,jj,kk)==0){
            voxelunion.grid(ii,jj,kk) = mm+1;
          }
        }
      }
    }
  }

  if(dim==3){
    std::ofstream out("vox_out.txt");
    if(!out.good()){
      std::cout<<"Cannot open output file\n";
      return -1;
    }
    out<<vx<<" "<<vy<<" "<<vz<<"\n";
    for(int ii = 0;ii<vx;ii++){
      for(int jj = 0;jj<vy;jj++){
        for(int kk = 0;kk<vz;kk++){
          out<<voxelunion.grid.read(ii,jj,kk)<<" ";
        }
        out<<"\n";
      }
      out<<"\n";
    }
    out.close();
    voxelunion.saveObj("voxel_out.obj");
  }

  //voxelunion.computeConnectivity();
  //if(dim==3){
    //voxelunion.saveConnectivity("fem.txt", true, false);
  //}else{
    //voxelunion.save2DSlice("quads.txt", slice, true);
  //}

  return 0;
}
