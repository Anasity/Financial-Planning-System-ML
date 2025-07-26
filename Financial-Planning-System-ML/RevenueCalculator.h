#pragma once // ������ �� ������������� ��������� ������������� �����

#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>
#include <msclr/marshal.h>

#define _USE_MATH_DEFINES
#include <math.h>

// ����� ��� ������� ������� � ��������
class RevenueCalculator {
public:
    // ��������� ��� �������� ������
    struct RevenueData {
        int tickets; // ���������� �������
        double revenue; // ����� �����
        RevenueData(int t, double r) : tickets(t), revenue(r) {}
    };

    struct CombinationResult {
        std::vector<std::function<double(double)>> e_funcs; // ������� ��������
        std::vector<std::function<double(double)>> r_funcs; // ������� �������
        std::string desc;
        std::vector<double> coefficients1;
        std::vector<double> coefficients2;
        double expense_offset;
        double revenue_offset;
        std::vector<std::pair<double, double>> intersection_points;
        std::string zone;
        double total_distance;
        size_t num_points;
    };

    // ������� ��������
    static double e1(double t) { return log(t + 1); }
    static double e2(double t) { return cos(t * 2 * M_PI / 8000000); }
    static double e3(double t) { return sin(t * 2 * M_PI / 12000000); }
    static double e4(double t) { return t; }

    // ������� �������
    static double r1(double t) { return sin(2 * M_PI * t / 15000000); }
    static double r2(double t) { return sin(2 * M_PI * t / 23000000); }
    static double r3(double t) { return (t * t); }
    static double r4(double t) { return sin(2 * M_PI * t / 50000000); }
    static double r5(double t) {
        return pow(t, 1.5) + 1e5 * t * sin(t / 1e7);
    }
    static double r6(double t) {
        return pow(t, 1.5) + 1e5 * t * sin(2 * M_PI * t / 50000000);
    }
    static double r7(double t) {
        return t;
    }
    static double r8(double t) {
        return pow(t, 1.2) * (1.0 + 0.3 * sin(t / 3e7));
    }
    static double r9(double t) {
        return pow(t, 3) + 1e3 * t * sin(t / 1e6);
    }
    static double r10(double t) { return (1 + 0.5 * sin(2 * M_PI * t / 45000000)); }


    static std::vector<double> solve_system(const std::vector<double>& t, // ���������� �������
        const std::vector<double>& xy_ratio, // ������
        const std::vector<std::function<double(double)>>& functions); // �������� ������� 

    static double calculate_distances(const std::pair<double, double>& user_point,
        const std::vector<std::pair<double, double>>& breakeven_points); // ������ ��� ������� ���������� �� ������ ������������ �� ��


    static CombinationResult* evaluate_combination(
        const std::vector<double>& t1, const std::vector<double>& xy_ratio1,
        const std::vector<double>& t2, const std::vector<double>& xy_ratio2,
        const std::pair<double, double>& user_point,
        const std::vector<std::function<double(double)>>& e_funcs,
        const std::vector<std::function<double(double)>>& r_funcs,
        const std::string& desc);

    //����� ��� ��������� ���� ��������� ���������� �������
    static std::vector<std::tuple<std::vector<std::function<double(double)>>,
        std::vector<std::function<double(double)>>,
        std::string>> generate_all_combinations();

    // ����� ��� ������ ����������� ����������
    // ����� ��� ������ ����������� ����������
    static CombinationResult* select_optimal_combination(
        const std::vector<CombinationResult*>& results,
        const std::string& user_zone,
        const std::pair<double, double>& user_point);

    //����� ��������� ������ 
    static void process_data(
        const std::pair<double, double>& user_point,
        void* richTextBoxResult,
        void* formHandle,
        const std::string& user_zone_str);

    // ����� ��� ��������� �������� ������ ��������� 
    static void draw_full_regression_chart(
        const std::vector<double>& t1, const std::vector<double>& xy_ratio1,
        const std::vector<double>& t2, const std::vector<double>& xy_ratio2,
        const std::vector<double>& full_coefficients1,
        const std::vector<double>& full_coefficients2,
        double full_expense_offset,
        double full_revenue_offset,
        void* formHandle);

    static const std::vector<std::function<double(double)>> all_e_functions; // ��� ������� ��������
    static const std::vector<std::function<double(double)>> all_r_functions; // ��� ������� �������
    static std::pair<double, double> calculate_gradient(
        const std::pair<double, double>& user_point,
        const std::vector<std::pair<double, double>>& breakeven_points,
        bool move_toward_points); // true - � ������ (������), false - �� ����� (�������)

    static std::string generate_recommendations(
        const std::pair<double, double>& gradient,
        const std::pair<double, double>& user_point);

    static void add_gradient_to_plot(
        std::ofstream& gp,
        const std::pair<double, double>& user_point,
        const std::pair<double, double>& gradient,
        const std::string& color,
        const std::vector<double>& t_values);

    static bool has_break_even_left(const std::vector<std::pair<double, double>>& points, double user_x) {
        for (const auto& point : points) {
            if (point.first < user_x) return true;
        }
        return false;
    }

    static bool has_break_even_right(const std::vector<std::pair<double, double>>& points, double user_x) {
        for (const auto& point : points) {
            if (point.first > user_x) return true;
        }
        return false;
    }

};