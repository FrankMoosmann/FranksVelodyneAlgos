#ifndef PNGIMAGEPROJECTOR_H_
#define PNGIMAGEPROJECTOR_H_

#include <LidarImageProjector.hpp>

/*!
 * \class PNGImageProjector
 * \brief The class implements projections from regular "sphere"-like images
 *
 *
 */
class PNGImageProjector : public LidarImageProjector
{
public:
  PNGImageProjector(std::string cfgFile);
  PNGImageProjector(double horizStartAngleRAD, double horizStopAngleRAD, double vertStartAngleRAD, double vertStopAngleRAD, unsigned int imgHSize, unsigned int imgVSize);
  virtual ~PNGImageProjector();

  virtual unsigned int getImgHorizSize() const;
  virtual unsigned int getImgVertSize() const;
  virtual double getImgHorizAngleRAD() const;
  virtual double getImgVertAngleRAD() const;

  virtual bool getImageIndexYP(double yawRAD, double pitchRAD, int &hi, int &vi) const;
  virtual void get3DCoordRel(int hi, int vi, double dist, double &x, double &y, double &z) const;
  virtual void getRotMatrix(int hi, int vi, matrixTools::DMatrix &mat) const;

private:
  class Spacing {
  public:
    virtual ~Spacing() {};
    virtual bool getImageIndex(double angleRAD, int &idx) = 0;
    virtual double getAngleRAD(int idx) = 0;
  };
  class RegularSpacing : public Spacing {
  private:
    double startAngleRAD; // is always in range [0,2pi)
    double rangeRAD; // positive value
    double rangeSign; // +/-1, direction of rangeRAD
    double pixSizeRAD; // can be positive or negative (same sign as rangeSign)
  public:
    RegularSpacing(double startAngleRAD, double stopAngleRAD, unsigned int nbPix); // +, +/-, +
    virtual ~RegularSpacing() {};
    virtual bool getImageIndex(double angleRAD, int &idx);
    virtual double getAngleRAD(int idx);
  };
  class IrregularSpacing : public Spacing {
  private:
    std::vector<double> vAnglesRAD;
    double avgPixSizeRAD; // can be positive or negative
    double pixSizeSign; // +/-1, sign of avgPixSizeRAD
  public:
    template <class InputIteratorT>
    IrregularSpacing(InputIteratorT angleRADBegin, InputIteratorT angleRADEnd, unsigned int nbPix) {
      while (angleRADBegin != angleRADEnd) { vAnglesRAD.push_back(*angleRADBegin); ++angleRADBegin; }
      avgPixSizeRAD = (vAnglesRAD.back()-vAnglesRAD.front())/(double)nbPix;
      pixSizeSign = (avgPixSizeRAD < 0.0) ? -1.0 : 1.0;
      std::cout << "  set up irregular spacing with " << vAnglesRAD.front() << ".." << vAnglesRAD.back() << ", pix " << nbPix << std::endl;
    };
    virtual ~IrregularSpacing() {};
    virtual bool getImageIndex(double angleRAD, int &idx);
    virtual double getAngleRAD(int idx);
  };
  Spacing *hSpacing, *vSpacing;
  unsigned int imgHSize;
  unsigned int imgVSize;
  double horizOpeningAngleRAD; // always positive
  double vertOpeningAngleRAD; // always positive

  double *xUnitSphere;
  double *yUnitSphere;
  double *zUnitSphere;
  double *rotMat;

  void setupSphere();
};


#endif // PNGIMAGEPROJECTOR_H_
