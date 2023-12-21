#include "basalt/rd_spline.h"
#include "basalt/so3_spline.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdbool.h>
#include <armadillo>
#include <cmath>

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

    std::string t;
    double x, y, z;
    double qx, qy, qz, qw;

    size_t nlines = 0;
    while (file >> t >> x >> y >> z >> qw >> qx >> qy >> qz)
    {
        nlines++;
        std::shared_ptr<Sophus::SE3d> pose_ptr(new Sophus::SE3d(Eigen::Quaterniond(qw, qx, qy, qz), Eigen::Vector3d(x, y, z)));
        poses.push_back(pose_ptr);
    }
    // std::cout << "Read " << nlines << " lines from file " << fn << std::endl;

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

    // std::cout << "99999999999" << std::endl;

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
    }

    // std::cout << "88888888888888" << std::endl;

    arma::cx_mat linear_v_x_freq = arma::fft(linear_v_x);
    arma::cx_mat linear_v_y_freq = arma::fft(linear_v_y);
    arma::cx_mat linear_v_z_freq = arma::fft(linear_v_z);

    // std::cout << "7777777777777" << std::endl;

    for (arma::uword i = linear_v_x_freq.n_elem / keep_freq; i < linear_v_x_freq.n_elem; i++) {
        linear_v_x_freq(i) = 0;
    }
    for (arma::uword i = linear_v_y_freq.n_elem / keep_freq; i < linear_v_y_freq.n_elem; i++) {
        linear_v_y_freq(i) = 0;
    }
    for (arma::uword i = linear_v_z_freq.n_elem / keep_freq; i < linear_v_z_freq.n_elem; i++) {
        linear_v_z_freq(i) = 0;
    }

    // std::cout << "6666666666" << std::endl;

    linear_v_x = arma::real(arma::ifft(linear_v_x_freq));
    linear_v_y = arma::real(arma::ifft(linear_v_y_freq));
    linear_v_z = arma::real(arma::ifft(linear_v_z_freq));

    // std::cout << "55555555555555" << std::endl;

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

void calc_fscore(
    std::vector<std::shared_ptr<Eigen::Vector3d>> & trans_deriv1,
    std::vector<std::shared_ptr<Eigen::Vector3d>> & rot_deriv1,
    std::vector<std::shared_ptr<Eigen::Vector3d>> & trans_deriv2,
    std::vector<std::shared_ptr<Eigen::Vector3d>> & rot_deriv2,
    double trans_threshold,
    double rot_threshold,
    double *fscore_trans,
    double *fscore_rot
) {
    size_t index = 0UL;
    long trans_num = 0;
    long rot_num = 0;

    while (1)
    {
        if (index < trans_deriv1.size() && index < rot_deriv1.size() && index < trans_deriv2.size() && index < rot_deriv2.size())
        {
            double x1 = (*(trans_deriv1[index]))(0);
            double y1 = (*(trans_deriv1[index]))(1);
            double z1 = (*(trans_deriv1[index]))(2);
            double rx1 = (*(rot_deriv1[index]))(0);
            double ry1 = (*(rot_deriv1[index]))(1);
            double rz1 = (*(rot_deriv1[index]))(2);
            double x2 = (*(trans_deriv2[index]))(0);
            double y2 = (*(trans_deriv2[index]))(1);
            double z2 = (*(trans_deriv2[index]))(2);
            double rx2 = (*(rot_deriv2[index]))(0);
            double ry2 = (*(rot_deriv2[index]))(1);
            double rz2 = (*(rot_deriv2[index]))(2);
            double trans_val = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
            double rot_val = sqrt((rx1-rx2)*(rx1-rx2) + (ry1-ry2)*(ry1-ry2) + (rz1-rz2)*(rz1-rz2));
            if (trans_val <= trans_threshold)
            {
                trans_num++;
            }
            if (rot_val <= rot_threshold)
            {
                rot_num++;
            }
            // std::cout << "rot_val " << rot_val << std::endl;
        } else {
            break;
        }
        index++;
    }

    double precision_trans = trans_num * 1.0 / trans_deriv1.size();
    double precision_rot = rot_num * 1.0 / rot_deriv1.size();
    double recall_trans = trans_num * 1.0 / trans_deriv2.size();
    double recall_rot = rot_num * 1.0 / rot_deriv2.size();

    *fscore_trans = 0.0;
    if (precision_trans + recall_trans > 0)
    {
        *fscore_trans = (2 * precision_trans * recall_trans) / (precision_trans + recall_trans);
    }
    *fscore_rot = 0.0;
    if (precision_rot + recall_rot > 0)
    {
        *fscore_rot = (2 * precision_rot * recall_rot) / (precision_rot + recall_rot);
    }
}

int main(int argc, char **argv)
{
    if (argc != 8)
    {
        std::cout << "Usage: " << argv[0] << " [/path/to/input/poses/file1] [/path/to/input/poses/file2] [argument controlling low pass filter] [sampling interval ns] [threshold start (inclusive)] [threshold end (inclusive)] [threshold interval]" << std::endl;
        return 0;
    }
    // else {
    //     std::cout << "Using " << argv[1] << " as the input poses file1" << std::endl;
    //     std::cout << "Using " << argv[2] << " as the input poses file2" << std::endl;
    //     std::cout << "Keeping first 1/" << argv[3] << " in the frequency domain" << std::endl;
    // }

    std::vector<std::shared_ptr<Sophus::SE3d>> poses1 = readPoses(std::string(argv[1]));

    std::vector<std::shared_ptr<Sophus::SE3d>> poses2 = readPoses(std::string(argv[2]));

    std::vector<std::shared_ptr<Eigen::Vector3d>> trans_deriv1;
    std::vector<std::shared_ptr<Eigen::Vector3d>> rot_deriv1;

    // std::cout << "1111111111" << std::endl;
    // estimate_v(200000000L, poses1, std::stoi(argv[3]), trans_deriv1, rot_deriv1);
    estimate_v(std::stol(argv[4]), poses1, std::stoi(argv[3]), trans_deriv1, rot_deriv1);

    std::vector<std::shared_ptr<Eigen::Vector3d>> trans_deriv2;
    std::vector<std::shared_ptr<Eigen::Vector3d>> rot_deriv2;

    // std::cout << "2222222222" << std::endl;
    // estimate_v(200000000L, poses2, std::stoi(argv[3]), trans_deriv2, rot_deriv2);
    estimate_v(std::stol(argv[4]), poses2, std::stoi(argv[3]), trans_deriv2, rot_deriv2);

    // std::cout << "3333333333" << std::endl;
    size_t size = rot_deriv1.size();

    // double fscore_trans, fscore_rot;
    // for (int i = 1; i <= 20; i++)
    // {
    //     calc_fscore(trans_deriv1, rot_deriv1, trans_deriv2, rot_deriv2, i * 0.005, 1000, &fscore_trans, &fscore_rot);
    //     std::cout << fscore_trans << " ";
    // }
    // std::cout << std::endl;

    double fscore_trans, fscore_rot;
    double threshold_start = std::stod(argv[5]);
    double threshold_end = std::stod(argv[6]);
    double threshold_interval = std::stod(argv[7]);
    double fscore_area_trans, fscore_area_rot;
    std::vector<double> fscore_transes;
    std::vector<double> fscore_rots;
    int num = 0;
    for (double threshold = threshold_start; threshold <= threshold_end + 1e-5; threshold += threshold_interval)
    {
        calc_fscore(trans_deriv1, rot_deriv1, trans_deriv2, rot_deriv2, threshold, threshold, &fscore_trans, &fscore_rot);
        double x_axis_len = exp(-10.0*(threshold-threshold_interval*0.5))-exp(-10.0*(threshold+threshold_interval*0.5));
        fscore_area_trans += fscore_trans*x_axis_len;
        fscore_area_rot += fscore_rot*x_axis_len;
        fscore_transes.push_back(fscore_trans);
        fscore_rots.push_back(fscore_rot);
        num++;
    }
    std::cout << fscore_area_trans << ", ";
    for (int i = 0; i < num; i++)
    {
        std::cout << fscore_transes[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << fscore_area_rot << ", ";
    for (int i = 0; i < num; i++)
    {
        std::cout << fscore_rots[i] << ", ";
    }
    std::cout << std::endl;

    // std::cout << "Translation robustness score is " << fscore_trans << std::endl;
    // std::cout << "Rotation robustness score is " << fscore_rot << std::endl;

    return 0;
}