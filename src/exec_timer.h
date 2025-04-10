//
// Created by jakob on 9/14/24.
//

#ifndef CPP_EXEC_TIMER_H
#define CPP_EXEC_TIMER_H

#define TIME_SECTIONS

#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>
#include <algorithm>

#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>
#include <algorithm>

class FunctionTimer {
private:
    struct FunctionData {
        std::string name;
        std::chrono::nanoseconds duration;

        FunctionData(const std::string &n) : name(n), duration(0) {}
    };

    struct SectionData {
        std::chrono::high_resolution_clock::time_point start_time;
        std::chrono::nanoseconds total_duration{0};
        int count = 0;
    };

    std::vector<FunctionData> functions;
    std::unordered_map<std::string, SectionData> sections;

public:
    template<typename Func, typename... Args>
    void register_and_run(const std::string &name, Func &&func, Args &&... args) {
        functions.emplace_back(name);
        auto &current_func = functions.back();

        std::cout << "Starting execution of: " << name << std::endl << std::flush;

        auto start_time = std::chrono::high_resolution_clock::now();
        std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
        auto end_time = std::chrono::high_resolution_clock::now();

        current_func.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    }

    void print_function_summary() const {
        if (functions.empty()) {
            std::cout << "No functions were timed." << std::endl;
            return;
        }

        size_t max_name_length = std::max_element(functions.begin(), functions.end(),
                                                  [](const FunctionData &a, const FunctionData &b) {
                                                      return a.name.length() < b.name.length();
                                                  }
        )->name.length();

        std::cout << std::endl << std::left << std::setw(max_name_length + 2) << "Function Name" << "| Execution Time"
                  << std::endl;
        std::cout << std::string(max_name_length + 2, '-') << "+" << std::string(20, '-') << std::endl;

        for (const auto &func: functions) {
            std::cout << std::left << std::setw(max_name_length + 2) << func.name << "| "
                      << format_duration(func.duration) << std::endl;
        }
    }

    void start_section(const std::string& name) {
#ifdef TIME_SECTIONS
        sections[name].start_time = std::chrono::high_resolution_clock::now();
#endif
    }

    void end_section(const std::string& name) {
#ifdef TIME_SECTIONS
        auto end_time = std::chrono::high_resolution_clock::now();
        auto& section = sections[name];
        section.total_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - section.start_time);
        section.count++;
#endif
    }

    void print_section_summary() const {
        if (sections.empty()) {
            std::cout << "No sections were timed." << std::endl;
            return;
        }

        size_t max_name_length = std::max_element(sections.begin(), sections.end(),
                                                  [](const auto &a, const auto &b) {
                                                      return a.first.length() < b.first.length();
                                                  }
        )->first.length();

        std::cout << std::endl << std::left << std::setw(max_name_length + 2) << "Section Name"
                  << "| Average Time    | Total Time      | Count" << std::endl;
        std::cout << std::string(max_name_length + 2, '-') << "+"
                  << std::string(17, '-') << "+"
                  << std::string(17, '-') << "+"
                  << std::string(7, '-') << std::endl;

        std::vector<std::pair<std::string, SectionData>> sorted_sections(sections.begin(), sections.end());
        std::sort(sorted_sections.begin(), sorted_sections.end(),
                  [](const auto &a, const auto &b) {
                      return a.second.total_duration > b.second.total_duration;
                  });

        for (const auto &[name, data]: sorted_sections) {
            std::chrono::nanoseconds avg_duration = data.count > 0 ?
                                                    data.total_duration / data.count : std::chrono::nanoseconds(0);

            std::cout << std::left << std::setw(max_name_length + 2) << name << "| "
                      << std::setw(15) << format_duration(avg_duration) << "| "
                      << std::setw(15) << format_duration(data.total_duration) << "| "
                      << std::setw(5) << data.count << std::endl;
        }
    }

private:
    static std::string format_duration(const std::chrono::nanoseconds &duration) {
        const auto ns = duration.count();

        if (ns >= 60'000'000'000) {  // More than a minute
            const auto minutes = ns / 60'000'000'000;
            const auto seconds = (ns % 60'000'000'000) / 1'000'000'000;
            return std::to_string(minutes) + "m " + std::to_string(seconds) + "s";
        } else if (ns >= 1'000'000'000) {  // More than a second
            const auto seconds = ns / 1'000'000'000;
            const auto milliseconds = (ns % 1'000'000'000) / 1'000'000;
            return std::to_string(seconds) + "." + std::to_string(milliseconds) + "s";
        } else if (ns >= 1'000'000) {  // More than a millisecond
            return std::to_string(ns / 1'000'000) + "ms";
        } else if (ns >= 1'000) {  // More than a microsecond
            return std::to_string(ns / 1'000) + "µs";
        } else {
            return std::to_string(ns) + "ns";
        }
    }
};

#endif //CPP_EXEC_TIMER_H