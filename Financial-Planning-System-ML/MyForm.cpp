#define NOMINMAX
#include "MyForm.h"
#include <windows.h>
using namespace System;
using namespace System::Windows::Forms;
using namespace msclr::interop;

// определения функций 
const std::vector<std::function<double(double)>> RevenueCalculator::all_e_functions = {
    &RevenueCalculator::e1,
    &RevenueCalculator::e2,
    &RevenueCalculator::e3,
};

const std::vector<std::function<double(double)>> RevenueCalculator::all_r_functions = {
    &RevenueCalculator::r1,
    &RevenueCalculator::r2,
    &RevenueCalculator::r3,
    &RevenueCalculator::r4,
    &RevenueCalculator::r5,
    &RevenueCalculator::r6,
    &RevenueCalculator::r7,
    &RevenueCalculator::r8,
    &RevenueCalculator::r9,
    &RevenueCalculator::r10,
};

int main(array<String^>^ args) {
    Application::EnableVisualStyles(); // визуальные стили
    Application::SetCompatibleTextRenderingDefault(false); // отображение текста
    CppCLRWinFormsProject::MyForm form; // создание формы приложения
    Application::Run(% form); // запуск приложения
    return 0;
}

// МЕТОД ДЛЯ РЕШЕНИЯ СИСТЕМЫ УРАВНЕНИЙ
std::vector<double> RevenueCalculator::solve_system(const std::vector<double>& t,
    const std::vector<double>& xy_ratio,
    const std::vector<std::function<double(double)>>& functions) {
    size_t N = t.size(); // количество месяцев
    size_t num_functions = functions.size(); // количество функций 
    // создание матрицы A и вектора b
    std::vector<std::vector<double>> A(num_functions, std::vector<double>(num_functions, 0.0));
    std::vector<double> b(num_functions, 0.0);

    for (size_t i = 0; i < num_functions; ++i) {
        for (size_t j = 0; j < num_functions; ++j) {
            double sum_val = 0.0;
            double alpha = 1e-10;
            for (size_t i = 0; i < num_functions; ++i) {
                A[i][i] += alpha;
            }

            for (size_t k = 0; k < N; ++k) {
                sum_val += functions[i](t[k]) * functions[j](t[k]);
            }
            A[i][j] = sum_val; // заполняем матрицу
        }
        double sum_val = 0.0;
        for (size_t k = 0; k < N; ++k) {
            sum_val += xy_ratio[k] * functions[i](t[k]);
        }
        b[i] = sum_val; // заполняем вектор b
    }
    // количество функций
    size_t n = num_functions;
    // находим коэффициенты
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i] / A[i][i]; //  находим значения 
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= A[k][i] * x[i]; // обновляем вектор b
        }
    }
    return x; // возвращаем найденные коэффициенты
}

// РАССЧЕТ РАССТОЯНИЙ ОТ ТОЧКИ ПОЛЬЗОВАТЕЛЯ ДО ТОЧЕК ПЕРЕСЕЧЕНИЯ
double RevenueCalculator::calculate_distances(
    const std::pair<double, double>& user_point,
    const std::vector<std::pair<double, double>>& breakeven_points) {
    if (breakeven_points.empty()) {
        return 0.0;
    }

    double total_distance = 0.0;
    for (const auto& point : breakeven_points) {
        double dx = user_point.first - point.first;
        double dy = user_point.second - point.second;
        total_distance += sqrt(dx * dx + dy * dy);
    }

    // Возвращаем среднее расстояние
    return total_distance / breakeven_points.size();
}

// МЕТОД ДЛЯ ОЦЕНКИ КОМБИНАЦИИ ФУНКЦИЙ
RevenueCalculator::CombinationResult* RevenueCalculator::evaluate_combination(
    const std::vector<double>& t1, const std::vector<double>& xy_ratio1,
    const std::vector<double>& t2, const std::vector<double>& xy_ratio2,
    const std::pair<double, double>& user_point,
    const std::vector<std::function<double(double)>>& e_funcs,
    const std::vector<std::function<double(double)>>& r_funcs,
    const std::string& desc) {

    try {
        // решение систем уравнений для расходов и доходов
        auto coefficients1 = solve_system(t1, xy_ratio1, e_funcs);
        auto coefficients2 = solve_system(t2, xy_ratio2, r_funcs);

        // считаем смещение для расходов
        double first_expense_point = xy_ratio1[0]; // первое реальное значение расходов
        double first_expense_t = t1[0]; // первое значение времени/количества
        double regression_value_at_first_point = 0.0;
        // считаем значение регрессии в первой точке
        for (size_t j = 0; j < e_funcs.size(); ++j) {
            // смещение = разница между реальным и предсказанным значением
            regression_value_at_first_point += coefficients1[j] * e_funcs[j](first_expense_t);
        }
        double expense_offset = first_expense_point - regression_value_at_first_point;

        // считаем смещение для доходов
        double first_revenue_point = xy_ratio2[0];
        double first_revenue_t = t2[0];
        regression_value_at_first_point = 0.0;
        for (size_t j = 0; j < r_funcs.size(); ++j) {
            regression_value_at_first_point += coefficients2[j] * r_funcs[j](first_revenue_t);
        }
        double revenue_offset = first_revenue_point - regression_value_at_first_point;


        // поиск точек пересечения с учетом смещений
        std::vector<double> t1_sorted = t1;
        std::sort(t1_sorted.begin(), t1_sorted.end()); // сортируем время

        std::vector<std::pair<double, double>> intersection_points;// вектор хранения тб

        for (size_t i = 1; i < t1_sorted.size(); ++i) {
            double t_prev = t1_sorted[i - 1];
            double t_curr = t1_sorted[i];

            double y_fitted1_prev = 0.0, y_fitted2_prev = 0.0;
            double y_fitted1_curr = 0.0, y_fitted2_curr = 0.0;

            for (size_t j = 0; j < e_funcs.size(); ++j) {
                y_fitted1_prev += coefficients1[j] * e_funcs[j](t_prev);
                y_fitted1_curr += coefficients1[j] * e_funcs[j](t_curr);
            }
            y_fitted1_prev += expense_offset;
            y_fitted1_curr += expense_offset;

            for (size_t j = 0; j < r_funcs.size(); ++j) {
                y_fitted2_prev += coefficients2[j] * r_funcs[j](t_prev);
                y_fitted2_curr += coefficients2[j] * r_funcs[j](t_curr);
            }
            y_fitted2_prev += revenue_offset;
            y_fitted2_curr += revenue_offset;

            // если функции пересекаются, вычисляем тб
            if ((y_fitted1_prev < y_fitted2_prev && y_fitted1_curr > y_fitted2_curr) ||
                (y_fitted1_prev > y_fitted2_prev && y_fitted1_curr < y_fitted2_curr)) {

                double t_intersect = t_prev + (t_curr - t_prev) * (y_fitted2_prev - y_fitted1_prev) /
                    ((y_fitted1_curr - y_fitted1_prev) - (y_fitted2_curr - y_fitted2_prev));

                double y_intersect = y_fitted1_prev + (y_fitted1_curr - y_fitted1_prev) *
                    (t_intersect - t_prev) / (t_curr - t_prev);

                intersection_points.emplace_back(t_intersect, y_intersect); // сохраняем тб
            }

        }

        if (intersection_points.empty()) {
            return nullptr;
        }

        // Рассчитываем общую дистанцию до точек пересечения
        double total_distance = calculate_distances(user_point, intersection_points);

        CombinationResult* result = new CombinationResult();
        result->e_funcs = e_funcs;
        result->r_funcs = r_funcs;
        result->desc = desc;
        result->coefficients1 = coefficients1;
        result->coefficients2 = coefficients2;
        result->expense_offset = expense_offset;
        result->revenue_offset = revenue_offset;
        result->intersection_points = intersection_points;
        result->total_distance = total_distance;
        result->num_points = intersection_points.size();

        // Расчет зоны для точки пользователя
        double user_t = user_point.first;
        double y1 = 0.0, y2 = 0.0;

        // Расчет значения расходов в точке пользователя
        for (size_t j = 0; j < e_funcs.size(); ++j) {
            y1 += coefficients1[j] * e_funcs[j](user_t);
        }
        y1 += expense_offset;

        // Расчет значения доходов в точке пользователя
        for (size_t j = 0; j < r_funcs.size(); ++j) {
            y2 += coefficients2[j] * r_funcs[j](user_t);
        }
        y2 += revenue_offset;

        // Определение зоны
        if (y2 > y1) {
            result->zone = "прибыль";
        }
        else {
            result->zone = "убыток";
        }

        return result;
    }
    catch (System::Exception^ e) {
        return nullptr;
    }
}

// Метод для генерации всех комбинаций функций
std::vector<std::tuple<std::vector<std::function<double(double)>>,
    std::vector<std::function<double(double)>>,
    std::string>> RevenueCalculator::generate_all_combinations() {

    const auto& all_e_functions = RevenueCalculator::all_e_functions;
    const auto& all_r_functions = RevenueCalculator::all_r_functions;

    std::vector<std::tuple<std::vector<std::function<double(double)>>,
        std::vector<std::function<double(double)>>,
        std::string>> combos; // вектор для храненния комбинаций

    // генерируем все подмножества для функций расходов
    for (size_t r = 1; r <= all_e_functions.size(); ++r) {
        std::vector<bool> mask(all_e_functions.size(), false);
        std::fill(mask.begin(), mask.begin() + r, true);

        do {
            std::vector<std::function<double(double)>> e_subset;
            std::string e_names;
            for (size_t i = 0; i < all_e_functions.size(); ++i) {
                if (mask[i]) {
                    e_subset.push_back(all_e_functions[i]);
                    if (!e_names.empty()) e_names += ", ";
                    e_names += "e" + std::to_string(i + 1); // название функции
                }
            }

            // генерируем все подмножества для функций доходов
            for (size_t s = 1; s <= all_r_functions.size(); ++s) {
                std::vector<bool> r_mask(all_r_functions.size(), false);
                std::fill(r_mask.begin(), r_mask.begin() + s, true);

                do {
                    std::vector<std::function<double(double)>> r_subset;
                    std::string r_names;
                    for (size_t j = 0; j < all_r_functions.size(); ++j) {
                        if (r_mask[j]) {
                            r_subset.push_back(all_r_functions[j]);
                            if (!r_names.empty()) r_names += ", ";
                            r_names += "r" + std::to_string(j + 1); // название функции 
                        }
                    }

                    std::string desc = "E: " + e_names + " | R: " + r_names;
                    combos.emplace_back(e_subset, r_subset, desc);
                } while (std::prev_permutation(r_mask.begin(), r_mask.end()));
            }
        } while (std::prev_permutation(mask.begin(), mask.end()));

    }

    return combos;
}

// ВЫБОР ОПТИМАЛЬНОЙ КОМБИНАЦИИ ИЗ РЕЗУЛЬТАТОВ


RevenueCalculator::CombinationResult* RevenueCalculator::select_optimal_combination(
    const std::vector<CombinationResult*>& results,
    const std::string& user_zone,
    const std::pair<double, double>& user_point)
{
    if (results.empty()) {
        return nullptr;
    }

    // Фильтрация результатов:
    std::vector<CombinationResult*> filtered_results;
    for (auto* result : results) {
        if (result->num_points >= 1) {
            // Проверка положения точки пользователя по вертикальной оси (Y)
            double user_t = user_point.first;
            double y_expense = 0.0;
            for (size_t j = 0; j < result->e_funcs.size(); ++j) {
                y_expense += result->coefficients1[j] * result->e_funcs[j](user_t);
            }
            y_expense += result->expense_offset;

            double y_revenue = 0.0;
            for (size_t j = 0; j < result->r_funcs.size(); ++j) {
                y_revenue += result->coefficients2[j] * result->r_funcs[j](user_t);
            }
            y_revenue += result->revenue_offset;

            bool condition_met = false;
            if (user_zone == "прибыль") {
                // Для зоны прибыли: доходы должны быть выше расходов в точке пользователя
                condition_met = (y_revenue > y_expense) &&
                    has_break_even_left(result->intersection_points, user_point.first);
            }
            else {
                // Для зоны убытка: доходы должны быть ниже расходов в точке пользователя
                condition_met = (y_revenue < y_expense) &&
                    has_break_even_right(result->intersection_points, user_point.first);
            }

            if (condition_met) {
                filtered_results.push_back(result);
            }
        }
    }

    if (filtered_results.empty()) {
        return nullptr;
    }

    // Выбор оптимальной комбинации
    auto comparator = [user_zone](const CombinationResult* a, const CombinationResult* b) {
        if (user_zone == "убыток") {
            // Для убытка выбираем минимальное расстояние (чтобы точка была ближе)
            return a->total_distance < b->total_distance;
        }
        else {
            // Для прибыли выбираем максимальное расстояние (чтобы точка была дальше)
            return a->total_distance > b->total_distance;
        }
        };

    auto it = std::min_element(filtered_results.begin(), filtered_results.end(), comparator);
    return *it;
}


// Реализация методов расчета градиента и рекомендаций
std::pair<double, double> RevenueCalculator::calculate_gradient(
    const std::pair<double, double>& user_point,
    const std::vector<std::pair<double, double>>& breakeven_points,
    bool move_toward_points)
{
    if (breakeven_points.empty()) {
        return { 0.0, 0.0 };
    }

    // Находим ближайшую правую точку безубыточности
    std::pair<double, double> right_point;
    bool has_right = false;
    double min_right_distance = std::numeric_limits<double>::max();

    for (const auto& point : breakeven_points) {
        if (point.first > user_point.first) {
            double dx = point.first - user_point.first;
            double dy = point.second - user_point.second;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < min_right_distance) {
                right_point = point;
                has_right = true;
                min_right_distance = dist;
            }
        }
    }

    double dQ = 0.0;
    double dP = 0.0;

    if (move_toward_points) {
        // Для зоны убытка - двигаемся только к ближайшей правой точке
        if (has_right) {
            double dx = right_point.first - user_point.first;
            double dy = right_point.second - user_point.second;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist > 0) {
                dQ += dx / dist;
                dP += dy / dist;
            }
        }
    }
    else {
        // Для зоны прибыли - двигаемся от всех точек безубыточности (старое поведение)
        // Находим ближайшие точки слева и справа
        std::pair<double, double> left_point, right_point;
        bool has_left = false;
        bool has_right = false;

        for (const auto& point : breakeven_points) {
            if (point.first < user_point.first) {
                if (!has_left || point.first > left_point.first) {
                    left_point = point;
                    has_left = true;
                }
            }
            else {
                if (!has_right || point.first < right_point.first) {
                    right_point = point;
                    has_right = true;
                }
            }
        }

        if (has_left) {
            double dx = user_point.first - left_point.first;
            double dy = user_point.second - left_point.second;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist > 0) {
                dQ += dx / dist;
                dP += dy / dist;
            }
        }
        if (has_right) {
            double dx = user_point.first - right_point.first;
            double dy = user_point.second - right_point.second;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist > 0) {
                dQ += dx / dist;
                dP += dy / dist;
            }
        }
    }

    // Нормализуем градиент
    double grad_magnitude = sqrt(dQ * dQ + dP * dP);
    if (grad_magnitude < 1e-10) {
        // Направление по умолчанию, если градиент нулевой
        dQ = (move_toward_points ? 1.0 : -1.0);
        dP = 0.0;
    }
    else {
        dQ /= grad_magnitude;
        dP /= grad_magnitude;
    }

    return { dQ, dP };
}

std::string RevenueCalculator::generate_recommendations(
    const std::pair<double, double>& gradient,
    const std::pair<double, double>& user_point)
{
    double dQ = gradient.first; // dΔ/dQ̃
    double dP = gradient.second; // dΔ/dP̃

    std::stringstream recommendations;
    recommendations << "\n";
    recommendations << "Рекомендации:\n\n";

    // 1. Рекомендация по изменению цены при изменении количества товара на 1
    if (dP != 0) {
        double price_change_gr = fabs(dQ / dP); // |(dΔ/dQ̃)/(dΔ/dQP̃)|
        recommendations << "1. На единицу товара: ";
        recommendations << "(" << price_change_gr << ", " << (dP > 0 ? "+1" : "-1") << ")\n";
        recommendations << "При изменении количества товара:\n";
        recommendations << "Если " << (dQ > 0 ? "увеличить" : "уменьшить")
            << " количество на 1 единицу, цену нужно "
            << (dP > 0 ? "увеличить" : "уменьшить")
            << " на " << price_change_gr << " руб.\n\n";
    }

    // 2. Рекомендация по изменению количества товара при изменении цены на 1 рубль
    if (dQ != 0) {
        double quantity_change_gr = fabs(dP / dQ); // |(dΔ/dP̃)/(dΔ/dQ̃)|
        recommendations << "2. На единицу цены: ";
        recommendations << "(" << quantity_change_gr << ", " << (dQ > 0 ? "+1" : "-1") << ")\n";
        recommendations << "При изменении цены:\n";
        recommendations << "Если " << (dP > 0 ? "увеличить" : "уменьшить")
            << " цену на 1 рубль, количество товара нужно "
            << (dQ > 0 ? "увеличить" : "уменьшить")
            << " на " << quantity_change_gr << " единиц.\n\n";
    }

    return recommendations.str();
}

void RevenueCalculator::add_gradient_to_plot(
    std::ofstream& gp,
    const std::pair<double, double>& user_point,
    const std::pair<double, double>& gradient,
    const std::string& color,
    const std::vector<double>& t_values)
{
    const double fixed_arrow_length = 6000000.0; // Фиксированная длина стрелки

    // Оцениваем масштаб осей (эмпирически или из данных)
    double x_range = 4e7; // Примерный диапазон по оси X (количество товара)
    double y_range = 1e11; // Примерный диапазон по оси Y (деньги)
    double scale_ratio = y_range / x_range; // Соотношение масштабов

    // Масштабируем вертикальную компоненту градиента
    double scaled_dP = gradient.second / scale_ratio;

    // Нормализуем градиент с учетом масштаба
    double grad_magnitude = sqrt(gradient.first * gradient.first + scaled_dP * scaled_dP);
    double end_x, end_y;

    if (grad_magnitude < 1e-10) {
        // Градиент нулевой - направление по умолчанию (например, вправо)
        end_x = user_point.first + fixed_arrow_length;
        end_y = user_point.second;
    }
    else {
        // Нормализованный градиент (с учетом масштаба)
        double normalized_dQ = gradient.first / grad_magnitude;
        double normalized_dP = scaled_dP / grad_magnitude;

        // Умножаем на фиксированную длину и возвращаем масштаб для Y
        end_x = user_point.first + normalized_dQ * fixed_arrow_length;
        end_y = user_point.second + normalized_dP * scale_ratio * fixed_arrow_length;
    }

    // Отрисовка стрелки
    gp << "set arrow from " << user_point.first << "," << user_point.second
        << " to " << end_x << "," << end_y
        << " head size 0.5,90 lw 2 lc rgb '" << color << "' front\n";
}


static System::Drawing::Bitmap^ GenerateGnuplotImage(
    const std::vector<double>& t1, const std::vector<double>& xy_ratio1,
    const std::vector<double>& t2, const std::vector<double>& xy_ratio2,
    const RevenueCalculator::CombinationResult* optimal_result,
    const std::pair<double, double>& user_point,
    const std::vector<double>& full_coefficients1,
    const std::vector<double>& full_coefficients2,
    double full_expense_offset,
    double full_revenue_offset,
    const std::pair<double, double>& gradient,
    const std::string& user_zone_str)
{
    std::string zone = (user_zone_str == "прибыль") ? "прибыль" : "убыток";
    // Удаление старых файлов
    std::remove("exp_data.txt");
    std::remove("rev_data.txt");
    std::remove("exp_fit.txt");
    std::remove("rev_fit.txt");
    std::remove("opt_exp_fit.txt");
    std::remove("opt_rev_fit.txt");
    std::remove("plot.gp");
    std::remove("chart.png");

    // Генерация уникального имени файла для графика
    std::string chart_filename = "chart_" + std::to_string(time(nullptr)) + ".png";

    // Запись данных в файлы (существующий код)
    std::ofstream exp_data("exp_data.txt");
    for (size_t i = 0; i < t1.size(); ++i) {
        exp_data << t1[i] << " " << xy_ratio1[i] << "\n";
    }
    exp_data.close();
    exp_data.close();

    std::ofstream rev_data("rev_data.txt");
    for (size_t i = 0; i < t2.size(); ++i) {
        rev_data << t2[i] << " " << xy_ratio2[i] << "\n";
    }
    rev_data.close();

    // Сортировка и запись регрессий
    std::vector<double> t1_sorted = t1;
    std::sort(t1_sorted.begin(), t1_sorted.end());

    std::ofstream exp_fit("exp_fit.txt");
    for (double t : t1_sorted) {
        double val = 0.0;
        for (size_t j = 0; j < RevenueCalculator::all_e_functions.size(); ++j) {
            val += full_coefficients1[j] * RevenueCalculator::all_e_functions[j](t);
        }
        val += full_expense_offset;
        exp_fit << t << " " << val << "\n";
    }
    exp_fit.close();

    std::ofstream rev_fit("rev_fit.txt");
    for (double t : t1_sorted) {
        double val = 0.0;
        for (size_t j = 0; j < RevenueCalculator::all_r_functions.size(); ++j) {
            val += full_coefficients2[j] * RevenueCalculator::all_r_functions[j](t);
        }
        val += full_revenue_offset;
        rev_fit << t << " " << val << "\n";
    }
    rev_fit.close();

    // Запись оптимальных регрессий, если есть
    if (optimal_result) {
        std::ofstream opt_exp_fit("opt_exp_fit.txt");
        for (double t : t1_sorted) {
            double val = 0.0;
            for (size_t j = 0; j < optimal_result->e_funcs.size(); ++j) {
                val += optimal_result->coefficients1[j] * optimal_result->e_funcs[j](t);
            }
            val += optimal_result->expense_offset;
            opt_exp_fit << t << " " << val << "\n";
        }
        opt_exp_fit.close();

        std::ofstream opt_rev_fit("opt_rev_fit.txt");
        for (double t : t1_sorted) {
            double val = 0.0;
            for (size_t j = 0; j < optimal_result->r_funcs.size(); ++j) {
                val += optimal_result->coefficients2[j] * optimal_result->r_funcs[j](t);
            }
            val += optimal_result->revenue_offset;
            opt_rev_fit << t << " " << val << "\n";
        }
        opt_rev_fit.close();

    }


    /// Создание скрипта gnuplot
    std::ofstream gp("plot.gp");
    gp << "set terminal pngcairo enhanced size 1200,700 font 'Arial,20'\n"
        << "set output '" << chart_filename << "'\n"

        << "set xlabel 'Количество товара (шт)'\n"
        << "set ylabel 'Деньги (руб)'\n"
        << "set key top left \n"
        << "set grid\n\n"
        << "# Вертикальная линия для точки пользователя\n"
        << "set arrow from " << user_point.first << ",graph(0) to "
        << user_point.first << ",graph(1) nohead lw 2 lc rgb '#AA00AA' dashtype 2\n\n";
    gp << "set rmargin 5\n"; // Увеличиваем правый отступ
    if (gradient.first != 0 || gradient.second != 0) {
        RevenueCalculator::add_gradient_to_plot(gp, user_point, gradient,
            user_zone_str == "убыток" ? "black" : "#black", t1);
    }

    // Если есть оптимальная комбинация и точки безубыточности
    if (optimal_result && !optimal_result->intersection_points.empty()) {
        // Сортируем точки безубыточности по X
        auto sorted_points = optimal_result->intersection_points;
        std::sort(sorted_points.begin(), sorted_points.end());

        // Добавляем закрашенные области между точками безубыточности
        // Сначала добавляем область до первой точки
        if (!t1_sorted.empty()) {
            double x0 = t1_sorted.front();
            double x1 = sorted_points.front().first;

            // Определяем цвет области перед первой точкой
            double mid_x = (x0 + x1) / 2;
            double mid_y_expense = 0.0;
            for (size_t j = 0; j < optimal_result->e_funcs.size(); ++j) {
                mid_y_expense += optimal_result->coefficients1[j] * optimal_result->e_funcs[j](mid_x);
            }
            mid_y_expense += optimal_result->expense_offset;

            double mid_y_revenue = 0.0;
            for (size_t j = 0; j < optimal_result->r_funcs.size(); ++j) {
                mid_y_revenue += optimal_result->coefficients2[j] * optimal_result->r_funcs[j](mid_x);
            }
            mid_y_revenue += optimal_result->revenue_offset;

            std::string zone_color = (mid_y_revenue > mid_y_expense) ? "green" : "red";

            gp << "set object " << 1 << " rect from " << x0 << ",graph(0) to "
                << x1 << ",graph(1) fc rgb '" << zone_color << "' fs solid 0.05 noborder\n";
        }

        // Затем добавляем области между точками
        for (size_t i = 0; i < sorted_points.size(); ++i) {
            double x1 = sorted_points[i].first;
            double y1 = sorted_points[i].second;

            // Определяем цвет области в зависимости от зоны справа от точки
            std::string zone_color = "red"; // Значение по умолчанию

            if (i < sorted_points.size() - 1) {
                // Для всех точек кроме последней
                double mid_x = (x1 + sorted_points[i + 1].first) / 2;

                double mid_y_expense = 0.0;
                for (size_t j = 0; j < optimal_result->e_funcs.size(); ++j) {
                    mid_y_expense += optimal_result->coefficients1[j] * optimal_result->e_funcs[j](mid_x);
                }
                mid_y_expense += optimal_result->expense_offset;

                double mid_y_revenue = 0.0;
                for (size_t j = 0; j < optimal_result->r_funcs.size(); ++j) {
                    mid_y_revenue += optimal_result->coefficients2[j] * optimal_result->r_funcs[j](mid_x);
                }
                mid_y_revenue += optimal_result->revenue_offset;

                zone_color = (mid_y_revenue > mid_y_expense) ? "green" : "red";
            }
            else {
                // Для последней точки (или единственной точки)
                double x2 = x1 + (t1_sorted.empty() ? x1 * 0.1 : (t1_sorted.back() - t1_sorted.front()) * 0.1);

                double y_expense = 0.0;
                for (size_t j = 0; j < optimal_result->e_funcs.size(); ++j) {
                    y_expense += optimal_result->coefficients1[j] * optimal_result->e_funcs[j](x2);
                }
                y_expense += optimal_result->expense_offset;

                double y_revenue = 0.0;
                for (size_t j = 0; j < optimal_result->r_funcs.size(); ++j) {
                    y_revenue += optimal_result->coefficients2[j] * optimal_result->r_funcs[j](x2);
                }
                y_revenue += optimal_result->revenue_offset;

                zone_color = (y_revenue > y_expense) ? "green" : "red";
            }

            // Добавляем закрашенную область
            if (!t1_sorted.empty()) {
                gp << "set object " << (i + 2) << " rect from " << x1 << ",graph(0) to "
                    << (i < sorted_points.size() - 1 ? sorted_points[i + 1].first : t1_sorted.back())
                    << ",graph(1) fc rgb '" << zone_color << "' fs solid 0.05 noborder\n";
            }
        }

        // Добавляем точки безубыточности с цветом в зависимости от зоны справа
        for (size_t i = 0; i < sorted_points.size(); ++i) {
            double x = sorted_points[i].first;
            double y = sorted_points[i].second;

            // Определяем цвет точки
            std::string point_color = "#FF0000"; // По умолчанию красный

            if (i < sorted_points.size() - 1) {
                double mid_x = (x + sorted_points[i + 1].first) / 2;

                double mid_y_expense = 0.0;
                for (size_t j = 0; j < optimal_result->e_funcs.size(); ++j) {
                    mid_y_expense += optimal_result->coefficients1[j] * optimal_result->e_funcs[j](mid_x);
                }
                mid_y_expense += optimal_result->expense_offset;

                double mid_y_revenue = 0.0;
                for (size_t j = 0; j < optimal_result->r_funcs.size(); ++j) {
                    mid_y_revenue += optimal_result->coefficients2[j] * optimal_result->r_funcs[j](mid_x);
                }
                mid_y_revenue += optimal_result->revenue_offset;

                point_color = (mid_y_revenue > mid_y_expense) ? "#00FF00" : "#FF0000";
            }
            else {
                // Для последней (или единственной) точки
                double x_right = x + (t1_sorted.empty() ? x * 0.1 : (t1_sorted.back() - t1_sorted.front()) * 0.1);

                double y_expense = 0.0;
                for (size_t j = 0; j < optimal_result->e_funcs.size(); ++j) {
                    y_expense += optimal_result->coefficients1[j] * optimal_result->e_funcs[j](x_right);
                }
                y_expense += optimal_result->expense_offset;

                double y_revenue = 0.0;
                for (size_t j = 0; j < optimal_result->r_funcs.size(); ++j) {
                    y_revenue += optimal_result->coefficients2[j] * optimal_result->r_funcs[j](x_right);
                }
                y_revenue += optimal_result->revenue_offset;

                point_color = (y_revenue > y_expense) ? "#00FF00" : "#FF0000";
            }

            gp << "set label " << (i + 1) << " at " << x << "," << y
                << " point pt 7 ps 2 lc rgb '" << point_color << "' front\n";
        }
    }

    gp << "plot ";

    // 1. тб (уже добавлены как labels)
    // 2. пользовательская точка
    gp << "'-' with points pt 7 ps 3 lc rgb '#AA00AA' title 'Пользовательская точка', \\\n";

    // 3. линии регрессий
    gp << "'exp_data.txt' with lines lw 1.5 lc rgb 'pink' title 'Данные расходов', \\\n"
        << "'rev_data.txt' with lines lw 1.5 lc rgb 'light-blue' title 'Данные доходов'";

    if (optimal_result) {
        gp << ", \\\n"
            << "'opt_exp_fit.txt' with lines lw 2 lc rgb 'red' title 'Оптимальная регрессия расходов', \\\n"
            << "'opt_rev_fit.txt' with lines lw 2 lc rgb 'blue' title 'Оптимальная регрессия доходов'";
    }

    gp << "\n";

    // Запись пользовательской точки
    gp << user_point.first << " " << user_point.second << "\ne\n";

    gp.close();

    // Запуск gnuplot
    system("gnuplot plot.gp");

    // Загрузка изображения с новым именем файла
    System::Drawing::Bitmap^ bmp = nullptr;
    try {
        bmp = gcnew System::Drawing::Bitmap(gcnew System::String(chart_filename.c_str()));
        std::remove(chart_filename.c_str()); // Удаляем временный файл
    }
    catch (System::Exception^ e) {
        // Обработка ошибки
        System::Drawing::Graphics^ g = System::Drawing::Graphics::FromImage(bmp);
        g->Clear(System::Drawing::Color::White);
        g->DrawString("Ошибка загрузки графика",
            gcnew System::Drawing::Font("Arial", 12),
            System::Drawing::Brushes::Black, 20, 20);
        delete g;
    }


    return bmp;
}

// РИСОВАНИЕ ГРАФИКОВ
static void draw_chart(
    const std::vector<double>& t1, const std::vector<double>& xy_ratio1,
    const std::vector<double>& t2, const std::vector<double>& xy_ratio2,
    const RevenueCalculator::CombinationResult* optimal_result,
    const std::pair<double, double>& user_point,
    const std::vector<double>& full_coefficients1,
    const std::vector<double>& full_coefficients2,
    double full_expense_offset,
    double full_revenue_offset,
    void* formHandle,
    const std::string& user_zone_str,
    const std::pair<double, double>& gradient)

{

    std::string zone = (user_zone_str == "прибыль") ? "прибыль" : "убыток";
    // Получаем форму
    System::Windows::Forms::Form^ form =
        static_cast<System::Windows::Forms::Form^>(
            System::Runtime::InteropServices::Marshal::GetObjectForIUnknown(
                System::IntPtr(formHandle)));


    CppCLRWinFormsProject::MyForm^ myForm = (CppCLRWinFormsProject::MyForm^)form;
    int panelWidth = myForm->panelChart->ClientSize.Width;
    int panelHeight = myForm->panelChart->ClientSize.Height;

    // Создаем изображение с графиком от gnuplot
    System::Drawing::Bitmap^ chartBmp = GenerateGnuplotImage(
        t1, xy_ratio1, t2, xy_ratio2, optimal_result, user_point,
        full_coefficients1, full_coefficients2,
        full_expense_offset, full_revenue_offset,
        gradient,
        user_zone_str);

    // Создаем итоговое изображение с легендой
    System::Drawing::Bitmap^ finalBmp = gcnew System::Drawing::Bitmap(panelWidth, panelHeight);
    System::Drawing::Graphics^ g = System::Drawing::Graphics::FromImage(finalBmp);
    g->Clear(System::Drawing::Color::White);
    g->SmoothingMode = System::Drawing::Drawing2D::SmoothingMode::AntiAlias;

    // Рисуем график
    g->DrawImage(chartBmp, 0, 0, panelWidth - 220, panelHeight - 30);

    // Обновляем график на форме
    myForm->UpdateChart(finalBmp);

}

//ГРАФИК С ПОЛНОЙ РЕГРЕССИЕЙ
void RevenueCalculator::draw_full_regression_chart(
    const std::vector<double>& t1, const std::vector<double>& xy_ratio1,
    const std::vector<double>& t2, const std::vector<double>& xy_ratio2,
    const std::vector<double>& full_coefficients1,
    const std::vector<double>& full_coefficients2,
    double full_expense_offset,
    double full_revenue_offset,
    void* formHandle)
{
    // Удаление старых файлов
    std::remove("exp_data.txt");
    std::remove("rev_data.txt");
    std::remove("exp_fit.txt");
    std::remove("rev_fit.txt");
    std::remove("plot.gp");

    // Запись данных расходов в файл
    std::ofstream exp_data("exp_data.txt");
    for (size_t i = 0; i < t1.size(); ++i) {
        exp_data << t1[i] << " " << xy_ratio1[i] << "\n";
    }
    exp_data.close();

    // Запись данных доходов в файл
    std::ofstream rev_data("rev_data.txt");
    for (size_t i = 0; i < t2.size(); ++i) {
        rev_data << t2[i] << " " << xy_ratio2[i] << "\n";
    }
    rev_data.close();

    // Сортировка данных для построения кривых
    std::vector<double> t1_sorted = t1;
    std::sort(t1_sorted.begin(), t1_sorted.end());

    // Запись регрессии расходов в файл
    std::ofstream exp_fit("exp_fit.txt");
    for (double t : t1_sorted) {
        double val = 0.0;
        for (size_t j = 0; j < RevenueCalculator::all_e_functions.size(); ++j) {
            val += full_coefficients1[j] * RevenueCalculator::all_e_functions[j](t);
        }
        val += full_expense_offset;
        exp_fit << t << " " << val << "\n";
    }
    exp_fit.close();

    std::ofstream rev_fit("rev_fit.txt");
    for (double t : t1_sorted) {
        double val = 0.0;
        for (size_t j = 0; j < all_r_functions.size(); ++j) {
            val += full_coefficients2[j] * all_r_functions[j](t);

        }
        val += full_revenue_offset;
        rev_fit << t << " " << val << "\n";
    }
    rev_fit.close();

    // Запись каждой функции доходов отдельно
    std::vector<std::string> r_func_files;
    for (size_t i = 0; i < RevenueCalculator::all_r_functions.size(); ++i) {
        std::string filename = "r_func_" + std::to_string(i + 1) + ".txt";
        std::ofstream r_func_file(filename);
        for (double t : t1_sorted) {
            double val = full_coefficients2[i] * RevenueCalculator::all_r_functions[i](t);

            r_func_file << t << " " << val << "\n";
        }
        r_func_file.close();
        r_func_files.push_back(filename);
    }


    // Запись каждой функции расходов отдельно
    std::vector<std::string> e_func_files;
    for (size_t i = 0; i < RevenueCalculator::all_e_functions.size(); ++i) {
        std::string filename = "r_func_" + std::to_string(i + 1) + ".txt";
        std::ofstream e_func_file(filename);
        for (double t : t1_sorted) {
            double val = full_coefficients2[i] * RevenueCalculator::all_e_functions[i](t);

            e_func_file << t << " " << val << "\n";
        }
        e_func_file.close();
        e_func_files.push_back(filename);
    }

    // Создание скрипта gnuplot

    std::ofstream gp("plot.gp");
    gp << "set terminal wxt size 1000,800 enhanced font 'Arial,14'\n\n"
        << "set title 'Аппроксимация данных'\n"
        << "set xlabel 'Количество товара (шт)'\n"
        << "set ylabel 'Деньги (руб)'\n"
        << "set key top left\n"
        << "set grid\n\n"
        << "set rmargin 5\n"  // <--- Перенесено перед plot
        << "plot 'exp_data.txt' with lines lw 1.5 lc rgb 'light-blue' title 'Данные расходов', \\\n"
        << "     'rev_data.txt' with lines lw 1.5 lc rgb 'pink' title 'Данные доходов', \\\n"
        << "     'exp_fit.txt' with lines lw 2 lc rgb 'blue' title 'Регрессия расходов', \\\n"
        << "     'rev_fit.txt' with lines lw 2 lc rgb 'red' title 'Регрессия доходов'\\\n";

    /*// Добавление каждой функции доходов в график
    std::vector<std::string> colors = { "red", "green", "blue", "purple", "orange",
                                     "brown", "magenta", "cyan", "yellow", "gray" };
    std::vector<std::string> titles = { "r1", "r2", "r3", "r4", "r5",
                                      "r6", "r7", "r8", "r9", "r10" };

    for (size_t i = 0; i < r_func_files.size(); ++i) {
        gp << ", \\\n     '" << r_func_files[i] << "' with lines lw 1.5 lc rgb '"
            << colors[i % colors.size()] << "' title '" << titles[i] << "'";
    }


    std::vector<std::string> colors1 = { "Pink", "Olive", "Gold" };
    std::vector<std::string> titles1 = { "e1", "e2", "e3" };

    for (size_t i = 0; i < e_func_files.size(); ++i) {
        gp << ", \\\n     '" << e_func_files[i] << "' with lines lw 1.5 lc rgb '"
            << colors[i % colors.size()] << "' title '" << titles[i] << "'";
    }*/

    gp << "\n";
    gp.close();

    // Запуск gnuplot
    system("gnuplot -persist plot.gp");
    // Удаление временных файлов
    for (const auto& file : e_func_files) {
        std::remove(file.c_str());
    }
}



// реализация process_data
void RevenueCalculator::process_data(
    const std::pair<double, double>& user_point,
    void* richTextBoxResult,
    void* formHandle,
    const std::string& user_zone_str)
{
    // Приведение void* к RichTextBox^
    System::Windows::Forms::RichTextBox^ resultBox =
        static_cast<System::Windows::Forms::RichTextBox^>(
            System::Runtime::InteropServices::Marshal::GetObjectForIUnknown(
                System::IntPtr(richTextBoxResult)));

    // Двойная очистка на всякий случай
    resultBox->Text = "";
    resultBox->Clear();
    Application::DoEvents(); // Даем время на обновление UI

    try {
        // чтение данных из файлов
        std::ifstream data_file("data.txt");
        if (!data_file) throw std::runtime_error("Could not open data.txt");

        std::vector<int> months, ticket_counts;
        for (int i = 0; i < 60; ++i) {
            int month, tickets;
            data_file >> month >> tickets;
            months.push_back(month);
            ticket_counts.push_back(tickets);
        }
        data_file.close();

        std::ifstream sebest_file("Sebest.txt");
        if (!sebest_file) throw std::runtime_error("Could not open Sebest.txt");

        std::vector<double> sebest_costs;
        for (int i = 0; i < 5; ++i) {
            double cost;
            sebest_file >> cost;
            sebest_costs.push_back(cost);
        }
        sebest_file.close();

        std::ifstream com_file("com_exp.txt");
        if (!com_file) throw std::runtime_error("Could not open com_exp.txt");

        std::vector<double> commercial_costs;
        for (int i = 0; i < 5; ++i) {
            double cost;
            com_file >> cost;
            commercial_costs.push_back(cost);
        }
        com_file.close();

        std::ifstream man_exp_file("man_exp.txt");
        if (!man_exp_file) throw std::runtime_error("Could not open man_exp.txt");

        std::vector<double> managerial_costs;
        for (int i = 0; i < 5; ++i) {
            double cost;
            man_exp_file >> cost;
            managerial_costs.push_back(cost);
        }
        man_exp_file.close();

        std::ifstream nalog_file("nalog.txt");
        if (!nalog_file) throw std::runtime_error("Could not open nalog.txt");

        std::vector<double> tax_costs;
        for (int i = 0; i < 5; ++i) {
            double cost;
            nalog_file >> cost;
            tax_costs.push_back(cost);
        }
        nalog_file.close();

        std::ifstream price_file("price.txt");
        if (!price_file) throw std::runtime_error("Could not open price.txt");

        std::vector<double> price_changes;
        for (int i = 0; i < 12; ++i) {
            double change;
            price_file >> change;
            price_changes.push_back(change);
        }
        price_file.close();

        std::ifstream revenue_file("revenue.txt");
        if (!revenue_file) throw std::runtime_error("Could not open revenue.txt");

        std::vector<double> revenue;
        for (int i = 0; i < 5; ++i) {
            double rev;
            revenue_file >> rev;
            revenue.push_back(rev);
        }
        revenue_file.close();

        // Расчет коэффициентов A и B для себестоимости
        double t_sum = 0.0, y_sum = 0.0, yt_sum = 0.0, t_square_sum = 0.0;
        for (int i = 0; i < 5; ++i) {
            t_sum += (i + 1);
            y_sum += sebest_costs[i];
            yt_sum += (i + 1) * sebest_costs[i];
            t_square_sum += pow(i + 1, 2);
        }

        double A_sebest = (yt_sum - (y_sum * t_sum) / 5) / (t_square_sum - (t_sum * t_sum) / 5);
        double B_sebest = (y_sum / 5) - A_sebest * (t_sum / 5);

        // Вычисляем коэффициенты для комм расх
        y_sum = 0.0;
        yt_sum = 0.0;
        for (int i = 0; i < 5; ++i) {
            y_sum += commercial_costs[i];
            yt_sum += (i + 1) * commercial_costs[i];
        }

        double A_com = (yt_sum - (y_sum * t_sum) / 5) / (t_square_sum - (t_sum * t_sum) / 5);
        double B_com = (y_sum / 5) - A_com * (t_sum / 5);

        // для управленческих расходов
        y_sum = 0.0;
        yt_sum = 0.0;
        for (int i = 0; i < 5; ++i) {
            y_sum += managerial_costs[i];
            yt_sum += (i + 1) * managerial_costs[i];
        }

        double A_man = (yt_sum - (y_sum * t_sum) / 5) / (t_square_sum - (t_sum * t_sum) / 5);
        double B_man = (y_sum / 5) - A_man * (t_sum / 5);

        // для налогов
        y_sum = 0.0;
        yt_sum = 0.0;
        for (int i = 0; i < 5; ++i) {
            y_sum += tax_costs[i];
            yt_sum += (i + 1) * tax_costs[i];
        }

        double A_tax = (yt_sum - (y_sum * t_sum) / 5) / (t_square_sum - (t_sum * t_sum) / 5);
        double B_tax = (y_sum / 5) - A_tax * (t_sum / 5);


        std::vector<double> variable_costs(60, 0.0);
        std::vector<double> monthly_costs(60, 0.0);
        std::vector<double> managerial_expenses(60, 0.0);
        std::vector<double> tax_expenses(60, 0.0);

        std::vector<double> total_tickets_per_year(5, 0.0);

        // сумма билетов по каждому году
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 12; ++j) {
                if (i * 12 + j < 60) {
                    total_tickets_per_year[i] += ticket_counts[i * 12 + j];
                }
            }
        }

        // расчет себестоимости по месяцам
        for (int i = 0; i < 60; ++i) {
            int year = i / 12;
            double year_cost = sebest_costs[year];
            variable_costs[i] = (year_cost - B_sebest) * (ticket_counts[i] / total_tickets_per_year[year]);
            monthly_costs[i] = (B_sebest / 12) + variable_costs[i];
        }

        // рассчет коммерческих расходов по месяцам
        std::vector<double> commercial_expenses(60, 0.0);
        for (int i = 0; i < 60; ++i) {
            int year = i / 12;
            double year_cost = commercial_costs[year];
            double variable_commercial_costs = (year_cost - B_com) * (ticket_counts[i] / total_tickets_per_year[year]);
            commercial_expenses[i] = (B_com / 12) + variable_commercial_costs;
        }

        // упр расходы по месяцам
        for (int i = 0; i < 60; ++i) {
            int year = i / 12;
            double year_cost = managerial_costs[year];
            double variable_managerial_costs = (year_cost - B_man) * (ticket_counts[i] / total_tickets_per_year[year]);
            managerial_expenses[i] = (B_man / 12) + variable_managerial_costs;
        }

        // налоги по месяцам
        for (int i = 0; i < 60; ++i) {
            int year = i / 12;
            double year_cost = tax_costs[year];
            double variable_tax_costs = (year_cost - B_tax) * (ticket_counts[i] / total_tickets_per_year[year]);
            tax_expenses[i] = (B_tax / 12) + variable_tax_costs;
        }

        // рассчет общих расходов 
        std::vector<double> total_expenses(60, 0.0);
        for (int i = 0; i < 60; ++i) {
            total_expenses[i] = monthly_costs[i] + commercial_expenses[i] +
                managerial_expenses[i] + tax_expenses[i];
        }

        // Подготовка данных для решения системы
        std::vector<std::pair<int, double>> expenses_by_tickets; // Вектор для хранения пар(количество билетов, расходы)
        for (int i = 0; i < 60; ++i) {
            expenses_by_tickets.emplace_back(ticket_counts[i], total_expenses[i]);
        }

        // сортировка по количеству билетов
        std::sort(expenses_by_tickets.begin(), expenses_by_tickets.end());

        // расчет доходов по месяцам
        std::vector<double> weighted_tickets(12, 0.0);
        std::vector<double> total_revenue(60, 0.0);

        for (int year = 0; year < 5; ++year) {
            double year_weighted_total = 0.0;
            for (int month = 0; month < 12; ++month) {
                int index = year * 12 + month;
                weighted_tickets[month] = ticket_counts[index] * (1 + price_changes[month]);
                year_weighted_total += weighted_tickets[month];
            }

            for (int month = 0; month < 12; ++month) {
                int index = year * 12 + month;
                total_revenue[index] = revenue[year] * (weighted_tickets[month] / year_weighted_total);
            }
        }

        // данные о доходах
        std::vector<RevenueData> revenue_data;
        for (int i = 0; i < 60; ++i) {
            revenue_data.emplace_back(ticket_counts[i], total_revenue[i]);
        }

        // сортировка расходов по количеству билетов
        std::sort(revenue_data.begin(), revenue_data.end(),
            [](const RevenueData& a, const RevenueData& b) {
                return a.tickets < b.tickets;
            });

        // Подготовка данных для решения системы уравнений
        std::vector<double> t1, xy_ratio1;
        for (const auto& pair : expenses_by_tickets) {
            t1.push_back(pair.first);
            xy_ratio1.push_back(pair.second);
        }

        std::vector<double> t2, xy_ratio2;
        for (const auto& data : revenue_data) {
            t2.push_back(data.tickets);
            xy_ratio2.push_back(data.revenue);
        }


        // генерация всех комбинаций функций 
        auto combinations = generate_all_combinations();
        String^ message = "Всего комбинаций: " + combinations.size().ToString() + "\n";
        resultBox->AppendText(message);

        // оценка всех комбинаций
        std::vector<RevenueCalculator::CombinationResult*> results;
        for (size_t i = 0; i < combinations.size(); ++i) {
            const auto& combo = combinations[i];
            auto result = evaluate_combination(t1, xy_ratio1, t2, xy_ratio2, user_point,
                std::get<0>(combo), std::get<1>(combo), std::get<2>(combo));

            if (result) {
                results.push_back(result);
            }

            if ((i + 1) % 10000000000 == 0) {
                message = "Пройдено " + (i + 1).ToString() + " комбинаций...\n";
                resultBox->AppendText(message);
            }
        }

        const auto& all_e_functions = RevenueCalculator::all_e_functions;
        const auto& all_r_functions = RevenueCalculator::all_r_functions;

        auto coefficients1 = solve_system(t1, xy_ratio1, all_e_functions);
        auto coefficients2 = solve_system(t2, xy_ratio2, all_r_functions);

        // Вычисляем смещения для полного набора функций
        double first_expense_point = xy_ratio1[0];
        double first_expense_t = t1[0];
        double regression_value_at_first_point = 0.0;
        for (size_t j = 0; j < all_e_functions.size(); ++j) {
            regression_value_at_first_point += coefficients1[j] * all_e_functions[j](first_expense_t);
        }
        double expense_offset = first_expense_point - regression_value_at_first_point;

        double first_revenue_point = xy_ratio2[0];
        double first_revenue_t = t2[0];
        regression_value_at_first_point = 0.0;
        for (size_t j = 0; j < all_r_functions.size(); ++j) {
            regression_value_at_first_point += coefficients2[j] * all_r_functions[j](first_revenue_t);
        }
        double revenue_offset = first_revenue_point - regression_value_at_first_point;

        // Рассчитываем зону с учетом смещений
        std::vector<double> t1_sorted = t1;
        std::sort(t1_sorted.begin(), t1_sorted.end());

        size_t idx = 0;
        while (idx < t1_sorted.size() && t1_sorted[idx] < user_point.first) {
            ++idx;
        }

        double y1, y2;
        if (idx == 0) {
            y1 = 0.0;
            for (size_t j = 0; j < all_e_functions.size(); ++j) {
                y1 += coefficients1[j] * all_e_functions[j](t1_sorted[0]);
            }
            y1 += expense_offset;

            y2 = 0.0;
            for (size_t j = 0; j < all_r_functions.size(); ++j) {
                y2 += coefficients2[j] * all_r_functions[j](t1_sorted[0]);
            }
            y2 += revenue_offset;
        }
        else if (idx == t1_sorted.size()) {
            y1 = 0.0;
            for (size_t j = 0; j < all_e_functions.size(); ++j) {
                y1 += coefficients1[j] * all_e_functions[j](t1_sorted.back());
            }
            y1 += expense_offset;

            y2 = 0.0;
            for (size_t j = 0; j < all_r_functions.size(); ++j) {
                y2 += coefficients2[j] * all_r_functions[j](t1_sorted.back());
            }
            y2 += revenue_offset;
        }
        else {
            double x0 = t1_sorted[idx - 1];
            double x1 = t1_sorted[idx];

            double y1_0 = 0.0, y1_1 = 0.0;
            for (size_t j = 0; j < all_e_functions.size(); ++j) {
                y1_0 += coefficients1[j] * all_e_functions[j](x0);
                y1_1 += coefficients1[j] * all_e_functions[j](x1);
            }
            y1_0 += expense_offset;


            y1_1 += expense_offset;
            y1 = y1_0 + (y1_1 - y1_0) * (user_point.first - x0) / (x1 - x0);

            double y2_0 = 0.0, y2_1 = 0.0;
            for (size_t j = 0; j < all_r_functions.size(); ++j) {
                y2_0 += coefficients2[j] * all_r_functions[j](x0);
                y2_1 += coefficients2[j] * all_r_functions[j](x1);
            }
            y2_0 += revenue_offset;
            y2_1 += revenue_offset;
            y2 = y2_0 + (y2_1 - y2_0) * (user_point.first - x0) / (x1 - x0);
        }

        // Создаем результат для полного набора функций
        RevenueCalculator::CombinationResult full_result;
        full_result.e_funcs = all_e_functions;
        full_result.r_funcs = all_r_functions;
        full_result.desc = "Все функции";
        full_result.coefficients1 = coefficients1;
        full_result.coefficients2 = coefficients2;
        full_result.expense_offset = expense_offset;
        full_result.revenue_offset = revenue_offset;
        full_result.intersection_points = {}; // Можно оставить пустым или рассчитать
        full_result.zone = "";
        full_result.total_distance = 0;
        full_result.num_points = 0;

        // Выбор оптимальной комбинации
        auto* optimal_result = select_optimal_combination(results, user_zone_str, user_point);

        if (!optimal_result) {
            resultBox->AppendText("не найдено комбинаций с т.б.\n");
            return;
        }

        message = "\nОптимальная комбинация:\n";
        message += "Функции: " + gcnew String(optimal_result->desc.c_str()) + "\n";
        message += "Количество Т.Б.: " + optimal_result->num_points.ToString() + "\n";
       // message += "Расстояние: " + optimal_result->total_distance.ToString() + "\n";
        resultBox->AppendText(message);

        // После выбора оптимальной комбинации
        if (optimal_result) {
            // Расчет градиента (true - к точкам для убытка, false - от точек для прибыли)
            bool move_toward = (user_zone_str == "убыток");
            auto gradient = calculate_gradient(user_point, optimal_result->intersection_points, move_toward);

            // Генерация рекомендаций
            std::string recommendations = generate_recommendations(gradient, user_point);
            resultBox->AppendText(gcnew System::String(recommendations.c_str()));

            // Обновление графика с градиентом
            draw_chart(t1, xy_ratio1, t2, xy_ratio2, optimal_result, user_point,
                coefficients1, coefficients2, expense_offset, revenue_offset,
                formHandle, user_zone_str, gradient);
        }

        for (auto* result : results) {
            delete result;
        }

        // Вызов для отображения графика с полными регрессиями и данными
        draw_full_regression_chart(t1, xy_ratio1, t2, xy_ratio2,
            coefficients1, coefficients2, expense_offset, revenue_offset, formHandle);

    }
    catch (System::Exception^ e) {
        resultBox->AppendText("Ошибка: " + e->Message + "\n");;
    }
}