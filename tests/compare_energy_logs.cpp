#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct LogRow {
    double pe = 0.0;
    double total_energy = 0.0;
    double temperature = 0.0;
};

struct LogMetrics {
    std::vector<LogRow> rows;
    double final_pe = 0.0;
    double final_total_energy = 0.0;
    double final_temperature = 0.0;
    double energy_drift = 0.0;
};

LogRow parse_row(const std::string& line) {
    std::istringstream input(line);
    double step = 0.0;
    double time = 0.0;
    double ke = 0.0;
    LogRow row;
    if (!(input >> step >> time >> row.pe >> ke >> row.total_energy >> row.temperature)) {
        throw std::runtime_error("Malformed energy-log row: " + line);
    }
    return row;
}

LogMetrics read_metrics(const std::string& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open energy log: " + path);
    }

    bool found_row = false;
    LogRow first;
    LogRow last;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        const LogRow row = parse_row(line);
        if (!found_row) {
            first = row;
            found_row = true;
        }
        last = row;
    }

    if (!found_row) {
        throw std::runtime_error("Energy log does not contain data rows: " + path);
    }

    input.clear();
    input.seekg(0, std::ios::beg);

    std::vector<LogRow> rows;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        rows.push_back(parse_row(line));
    }

    return LogMetrics{
        .rows = std::move(rows),
        .final_pe = last.pe,
        .final_total_energy = last.total_energy,
        .final_temperature = last.temperature,
        .energy_drift = last.total_energy - first.total_energy,
    };
}

double parse_tolerance(const char* value, const char* name) {
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0' || parsed < 0.0) {
        throw std::runtime_error(std::string("Invalid ") + name + " tolerance: " + value);
    }
    return parsed;
}

void require_close(double serial,
                   double parallel,
                   double tolerance,
                   const char* metric) {
    const double error = std::abs(serial - parallel);
    if (error > tolerance) {
        std::ostringstream message;
        message << metric << " differs by " << error
                << " (serial=" << serial
                << ", mpi=" << parallel
                << ", tolerance=" << tolerance << ')';
        throw std::runtime_error(message.str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 6) {
            throw std::runtime_error(
                "Usage: gmd_compare_energy_logs serial.log mpi.log pe_tol temp_tol drift_tol");
        }

        const LogMetrics serial = read_metrics(argv[1]);
        const LogMetrics parallel = read_metrics(argv[2]);
        if (serial.rows.size() != parallel.rows.size()) {
            throw std::runtime_error(
                "Energy log row count differs between serial and MPI runs");
        }

        const double row_tolerance = parse_tolerance(argv[3], "per-step energy");
        const double temp_tolerance = parse_tolerance(argv[4], "final temperature");
        for (std::size_t index = 0; index < serial.rows.size(); ++index) {
            require_close(serial.rows[index].pe,
                          parallel.rows[index].pe,
                          row_tolerance,
                          "Per-step PE");
            require_close(serial.rows[index].total_energy,
                          parallel.rows[index].total_energy,
                          row_tolerance,
                          "Per-step total energy");
            require_close(serial.rows[index].temperature,
                          parallel.rows[index].temperature,
                          temp_tolerance,
                          "Per-step temperature");
        }

        require_close(serial.final_pe,
                      parallel.final_pe,
                      row_tolerance,
                      "Final PE");
        require_close(serial.final_total_energy,
                      parallel.final_total_energy,
                      row_tolerance,
                      "Final total energy");
        require_close(serial.final_temperature,
                      parallel.final_temperature,
                      temp_tolerance,
                      "Final temperature");
        require_close(serial.energy_drift,
                      parallel.energy_drift,
                      parse_tolerance(argv[5], "NVE energy drift"),
                      "NVE energy drift");
    } catch (const std::exception& error) {
        std::cerr << "MPI consistency check failed: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
