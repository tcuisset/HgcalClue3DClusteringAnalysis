#ifndef CLUSTER_H
#define CLUSTER_H

class Cluster {
  public:
    Cluster():energy_(0.f), x_(0.f), y_(0.f), z_(0.f), layer_(-1) {};
    ~Cluster() = default;

  void addCell(int idx) { hits_.push_back(idx);}

  void addEnergy(float value) {
    energy_ += value;
  }

  void addEnergyAndRescale(float value, float z) {
    // Do hard-coded rescaling...
    if (z < 54.f) {
      energy_ += value * 0.0105;
    } else if (z > 54.f and z < 154.f) {
      energy_ += value * 0.0812;
    } else if (z > 154.f) {
      energy_ += value * 0.12508;
    }
  }

  void setLayer(int value) {
    if (layer_ != -1) {
      assert(value==layer_);
    }
    layer_ = value;
  }

  void setPosition(float x, float y, float z) {
    x_ = x;
    y_ = y;
    z_ = z;
  }

  std::tuple<float, float, float> position() const {
    return std::make_tuple(x_, y_, z_);
  }

  float energy () const {return energy_;}
  int layer() const {return layer_;}
  const std::vector<int> & hits() const { return hits_;}

  private:
    std::vector<int> hits_;
    float energy_;
    float x_;
    float y_;
    float z_;

    int layer_;
};

#endif
