#include "basalt/rd_spline.h"
#include "basalt/so3_spline.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdbool.h>
#include <armadillo>

std::vector<std::shared_ptr<Sophus::SE3d>> readPoses(std::string fn)
{
    std::ifstream file;
    file.open(fn);

    std::vector<std::shared_ptr<Sophus::SE3d>> poses;

    if (!file.is_open())
    {
        std::cout << "Can't open file " << fn << std::endl;
        return poses;
    }

    unsigned secs, nsecs;
    double x, y, z;
    double qx, qy, qz, qw;
    char ch;

    size_t nlines = 0;
    while (file >> secs >> ch >> nsecs >> ch >> x >> ch >> y >> ch >> z >> ch >> qw >> ch >> qx >> ch >> qy >> ch >> qz)
    {
        nlines++;
        std::shared_ptr<Sophus::SE3d> pose_ptr(new Sophus::SE3d(Eigen::Quaterniond(qw, qx, qy, qz), Eigen::Vector3d(x, y, z)));
        poses.push_back(pose_ptr);
    }
    std::cout << "Read " << nlines << " lines from file " << fn << std::endl;

    file.close();

    return poses;
}

void estimate_v(long interv_ns, std::vector<std::shared_ptr<Sophus::SE3d>> & poses, int keep_freq, std::vector<std::shared_ptr<Eigen::Vector3d>> & trans_deriv, std::vector<std::shared_ptr<Eigen::Vector3d>> & rot_deriv)
{
    basalt::So3Spline<4> rot_spline(interv_ns);
    basalt::RdSpline<3, 4> pos_spline(interv_ns);

    for (std::shared_ptr<Sophus::SE3d> pose_ptr : poses) {
        rot_spline.knotsPushBack(pose_ptr->so3());
        pos_spline.knotsPushBack(pose_ptr->translation());
    }

    std::vector<std::shared_ptr<Sophus::SE3d>> poses_deriv;
    int n_poses = poses.size();

    arma::mat linear_v_x(n_poses - 6, 1);
    arma::mat linear_v_y(n_poses - 6, 1);
    arma::mat linear_v_z(n_poses - 6, 1);
    arma::mat angular_v_x(n_poses - 6, 1);
    arma::mat angular_v_y(n_poses - 6, 1);
    arma::mat angular_v_z(n_poses - 6, 1);

    for (size_t i = 3; i < n_poses - 3; i++)
    {
        Eigen::Vector3d linear_v = pos_spline.velocity(interv_ns * i);
        linear_v_x(i - 3) = linear_v(0);
        linear_v_y(i - 3) = linear_v(1);
        linear_v_z(i - 3) = linear_v(2);
        Eigen::Vector3d angular_v = rot_spline.velocityBody(interv_ns * i);
        angular_v_x(i - 3) = angular_v(0);
        angular_v_y(i - 3) = angular_v(1);
        angular_v_z(i - 3) = angular_v(2);

        // poses_deriv.push_back(std::shared_ptr<Sophus::SE3d>(new Sophus::SE3d(
        //     Sophus::SO3d::exp(rot_spline.velocityBody(interv_ns * i)),
        //     pos_spline.velocity(interv_ns * i)
        // )));
    }

    // std::cout << "Time domain v_x rows " << linear_v_x.n_rows << std::endl;
    // std::cout << "Time domain v_x cols " << linear_v_x.n_cols << std::endl;
    // std::cout << "Time domain v_y rows " << linear_v_y.n_rows << std::endl;
    // std::cout << "Time domain v_y cols " << linear_v_y.n_cols << std::endl;
    // std::cout << "Time domain v_z rows " << linear_v_z.n_rows << std::endl;
    // std::cout << "Time domain v_z cols " << linear_v_z.n_cols << std::endl;

    arma::cx_mat linear_v_x_freq = arma::fft(linear_v_x);
    arma::cx_mat linear_v_y_freq = arma::fft(linear_v_y);
    arma::cx_mat linear_v_z_freq = arma::fft(linear_v_z);

    // std::cout << "Frequency domain v_x rows " << linear_v_x_freq.n_rows << std::endl;
    // std::cout << "Frequency domain v_x cols " << linear_v_x_freq.n_cols << std::endl;
    // std::cout << "Frequency domain v_y rows " << linear_v_y_freq.n_rows << std::endl;
    // std::cout << "Frequency domain v_y cols " << linear_v_y_freq.n_cols << std::endl;
    // std::cout << "Frequency domain v_z rows " << linear_v_z_freq.n_rows << std::endl;
    // std::cout << "Frequency domain v_z cols " << linear_v_z_freq.n_cols << std::endl;

    // std::ofstream fvx, fvy, fvz;
    // fvx.open("../result/fvx.txt");
    // fvy.open("../result/fvy.txt");
    // fvz.open("../result/fvz.txt");

    // for (arma::uword i = 0; i < linear_v_x_freq.n_elem; i++) {
    //     fvx << std::real(linear_v_x_freq(i)) << "+j" << std::imag(linear_v_x_freq(i)) << std::endl;
    // }
    // for (arma::uword i = 0; i < linear_v_y_freq.n_elem; i++) {
    //     fvy << std::real(linear_v_y_freq(i)) << "+j" << std::imag(linear_v_y_freq(i)) << std::endl;
    // }
    // for (arma::uword i = 0; i < linear_v_z_freq.n_elem; i++) {
    //     fvz << std::real(linear_v_z_freq(i)) << "+j" << std::imag(linear_v_z_freq(i)) << std::endl;
    // }

    // fvx.close();
    // fvy.close();
    // fvz.close();

    for (arma::uword i = linear_v_x_freq.n_elem / keep_freq; i < linear_v_x_freq.n_elem; i++) {
        linear_v_x_freq(i) = 0;
    }
    for (arma::uword i = linear_v_y_freq.n_elem / keep_freq; i < linear_v_y_freq.n_elem; i++) {
        linear_v_y_freq(i) = 0;
    }
    for (arma::uword i = linear_v_z_freq.n_elem / keep_freq; i < linear_v_z_freq.n_elem; i++) {
        linear_v_z_freq(i) = 0;
    }

    linear_v_x = arma::real(arma::ifft(linear_v_x_freq));
    linear_v_y = arma::real(arma::ifft(linear_v_y_freq));
    linear_v_z = arma::real(arma::ifft(linear_v_z_freq));

    arma::cx_mat angular_v_x_freq = arma::fft(angular_v_x);
    arma::cx_mat angular_v_y_freq = arma::fft(angular_v_y);
    arma::cx_mat angular_v_z_freq = arma::fft(angular_v_z);

    for (arma::uword i = angular_v_x_freq.n_elem / keep_freq; i < angular_v_x_freq.n_elem; i++) {
        angular_v_x_freq(i) = 0;
    }
    for (arma::uword i = angular_v_y_freq.n_elem / keep_freq; i < angular_v_y_freq.n_elem; i++) {
        angular_v_y_freq(i) = 0;
    }
    for (arma::uword i = angular_v_z_freq.n_elem / keep_freq; i < angular_v_z_freq.n_elem; i++) {
        angular_v_z_freq(i) = 0;
    }

    angular_v_x = arma::real(arma::ifft(angular_v_x_freq));
    angular_v_y = arma::real(arma::ifft(angular_v_y_freq));
    angular_v_z = arma::real(arma::ifft(angular_v_z_freq));

    for (size_t i = 0; i < n_poses - 6; i++) {
        rot_deriv.push_back(std::shared_ptr<Eigen::Vector3d>(new Eigen::Vector3d(angular_v_x(i), angular_v_y(i), angular_v_z(i))));
        trans_deriv.push_back(std::shared_ptr<Eigen::Vector3d>(new Eigen::Vector3d(linear_v_x(i), linear_v_y(i), linear_v_z(i))));
    }
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cout << "Usage: " << argv[0] << " [/path/to/input/poses/file1] [/path/to/input/poses/file2] [argument controlling low pass filter]" << std::endl;
        return 0;
    } else {
        std::cout << "Using " << argv[1] << " as the input poses file1" << std::endl;
        std::cout << "Using " << argv[2] << " as the input poses file2" << std::endl;
        std::cout << "Keeping first 1/" << argv[3] << " in the frequency domain" << std::endl;
    }

    std::vector<std::shared_ptr<Sophus::SE3d>> poses1 = readPoses(std::string(argv[1]));

    std::vector<std::shared_ptr<Sophus::SE3d>> poses2 = readPoses(std::string(argv[2]));

    std::vector<std::shared_ptr<Eigen::Vector3d>> trans_deriv1;
    std::vector<std::shared_ptr<Eigen::Vector3d>> rot_deriv1;

    estimate_v(200000000L, poses, std::stoi(argv[3]), trans_deriv1, rot_deriv1);

    std::vector<std::shared_ptr<Eigen::Vector3d>> trans_deriv2;
    std::vector<std::shared_ptr<Eigen::Vector3d>> rot_deriv2;

    estimate_v(200000000L, poses, std::stoi(argv[3]), trans_deriv2, rot_deriv2);

    // TODO: calculate the robustness score based on the difference between (trans_deriv1, rot_deriv1) and (trans_deriv2, rot_deriv2)

    return 0;
}