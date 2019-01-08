#include <iostream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <LidarImageProjectorPNG.hpp>
#include <PngDistanceImage.hpp>

struct Point {
  float x;
  float y;
  float z;
  float reflectivity;
  Point(float x, float y, float z, float r) : x(x), y(y), z(z), reflectivity(r) {
  }
};

using Pointcloud = std::vector<Point>;

Pointcloud loadIntoPointcloud(std::string kittiVeloPointcloudFile) {
  Pointcloud pointcloud;
  // the following code on how to parse a kitti velodyne raw file is directly taken from the kitti documentation
  // allocate 4 MB buffer (only ~130*4*4 KB are needed)
  int32_t num = 1000000;
  float *data = (float*)malloc(num*sizeof(float));

  // pointers
  float *px = data+0;
  float *py = data+1;
  float *pz = data+2;
  float *pr = data+3;

  // load point cloud
  FILE *stream;
  stream = fopen (kittiVeloPointcloudFile.c_str(),"rb");
  num = fread(data,sizeof(float),num,stream)/4;
  for (int32_t i=0; i<num; i++) {
    pointcloud.push_back(Point(*px,*py,*pz,*pr));
    px+=4; py+=4; pz+=4; pr+=4;
  }
  fclose(stream);

  return pointcloud;
}

int main(int argc, char *argv[])
{
  using namespace std;
  namespace fs = boost::filesystem;
  namespace po = boost::program_options;

  string infile = "";
  string outfile = "";
  string imageConfigFile = "";

  po::options_description progOpt;
  progOpt.add_options()
    ("help", "produce help message")
    ("infile", po::value(&infile), "intput file: kitti binary velodyne pointcloud")
    ("outfile", po::value(&outfile), "output file: 16bit png distance image")
    ("configFile", po::value(&imageConfigFile), "file containing the configuration of the image projection")
    ;
  po::positional_options_description posOpt;
  posOpt.add("infile", 1);
  posOpt.add("outfile", 1);
  posOpt.add("configFile", 1);

  po::variables_map poVM;
  try {
    //po::store(parse_command_line(argc, argv, mainParams, po::command_line_style::unix_style ^ po::command_line_style::allow_short), poVM);
    po::store(po::command_line_parser(argc, argv).options(progOpt).positional(posOpt).run(), poVM);
    po::notify(poVM);
  } catch (exception &e) {
    cerr << "Error when parsing arguments:" << e.what() << endl;
    return 1;
  }

  if (poVM.count("help")) {
    cout << progOpt << endl;
    return 0;
  }

  if (!fs::exists(infile) || !fs::is_regular_file(infile)) {
    cout << "Input file \"" << infile << "\" does not exist or is not a file" << endl;
    return 1;
  }

  if (!fs::exists(imageConfigFile) || !fs::is_regular_file(imageConfigFile)) {
    cout << "Image config file \"" << infile << "\" does not exist or is not a file" << endl;
    return 1;
  }

  auto pointcloud = loadIntoPointcloud(infile);
  cout << "loaded " << pointcloud.size() << " points" << endl;

  PNGImageProjector projector(imageConfigFile);
  PngDistanceImage pngImage(projector.getImgHorizSize(), projector.getImgVertSize());
  pngImage.fill(0.0);
  for (auto &point : pointcloud) {
    int hi, vi;
    double dist;
    if (projector.getImageIndexRel(point.x, point.y, point.z, hi, vi, dist)) {
      pngImage.setDistance(vi, hi, dist);
    }
  }
  pngImage.save(outfile);
}
