#include "gmd/io/frame_writer.hpp"

#include <array>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gmd/core/physical_constants.hpp"
#include "gmd/integrator/thermostat.hpp"  // temperature_from_twice_ke
#include "gmd/system/system.hpp"

namespace gmd {

namespace {

constexpr std::array<const char*, 119> kElementSymbols = {
    "X",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
};

void append_real(std::string& out, double value) {
    char buffer[40];
    const int written = std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    out.append(buffer, static_cast<std::size_t>(written));
    // Keep every real recognisably real ("20" -> "20.0"), so extxyz readers
    // type a value the same way in every frame.
    bool has_marker = false;
    for (int i = 0; i < written; ++i) {
        const char c = buffer[i];
        if (c == '.' || c == 'e' || c == 'n' || c == 'i') {  // also nan / inf
            has_marker = true;
            break;
        }
    }
    if (!has_marker) out += ".0";
}

void append_key_real(std::string& out, const char* key, double value) {
    out.push_back(' ');
    out += key;
    out.push_back('=');
    append_real(out, value);
}

}  // namespace

const char* element_symbol(int atomic_number) noexcept {
    if (atomic_number <= 0 || atomic_number >= static_cast<int>(kElementSymbols.size())) {
        return kElementSymbols[0];
    }
    return kElementSymbols[static_cast<std::size_t>(atomic_number)];
}

FrameWriter::~FrameWriter() {
    close();
}

void FrameWriter::open(const std::filesystem::path& path, bool write_velocities,
                       FrameProvenance provenance) {
    close();
    out_.open(path, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!out_.is_open()) {
        throw std::runtime_error("FrameWriter: failed to open " + path.string());
    }
    write_velocities_ = write_velocities;
    provenance_ = std::move(provenance);
    frame_count_ = 0;
}

void FrameWriter::close() {
    if (out_.is_open()) {
        out_.flush();
        out_.close();
    }
}

std::string FrameWriter::format_frame(const System& system, std::uint64_t step, double time_fs,
                                      double twice_ke, std::size_t dof, bool write_velocities,
                                      const FrameProvenance& provenance) {
    const std::size_t n = system.atom_count();
    const auto& lengths = system.box().lengths;
    const double potential_energy = system.potential_energy();
    const double kinetic_energy = 0.5 * twice_ke;
    const double temperature = dof > 0 ? temperature_from_twice_ke(twice_ke, dof) : 0.0;
    const double volume = lengths[0] * lengths[1] * lengths[2];

    std::string out;
    out.reserve(512 + n * (write_velocities ? 240 : 170));
    out += std::to_string(n);
    out.push_back('\n');

    out += "Lattice=\"";
    append_real(out, lengths[0]);
    out += " 0 0 0 ";
    append_real(out, lengths[1]);
    out += " 0 0 0 ";
    append_real(out, lengths[2]);
    out += "\" Properties=id:I:1:species:S:1:pos:R:3:forces:R:3";
    if (write_velocities) out += ":vel:R:3";
    out += " pbc=\"T T T\"";
    out += " step=" + std::to_string(step);
    append_key_real(out, "time_fs", time_fs);
    append_key_real(out, "energy", potential_energy);
    append_key_real(out, "potential_energy", potential_energy);
    append_key_real(out, "kinetic_energy", kinetic_energy);
    append_key_real(out, "temperature", temperature);
    append_key_real(out, "volume", volume);
    out += " format=";
    out += kFrameFormat;
    out += " units=\"energy:eV length:Angstrom force:eV/Angstrom velocity:Angstrom/fs "
           "time:fs temperature:K\"";
    out += " run_id=" + (provenance.run_id.empty() ? std::string("none") : provenance.run_id);
    out += " force_field=" +
           (provenance.force_field_kind.empty() ? std::string("none")
                                                : provenance.force_field_kind);
    out += " model_sha256=" +
           (provenance.model_sha256.empty() ? std::string("none") : provenance.model_sha256);
    if (!provenance.model_changed_from.empty()) {
        out += " model_changed_from=" + provenance.model_changed_from;
    }
    out.push_back('\n');

    // Atom lines in TAG order. Serial runs never reorder atoms (tag == index),
    // and the MPI output system is assembled by tag, but sorting here keeps
    // the guarantee independent of either.
    std::vector<std::size_t> order(n);
    std::vector<std::size_t> by_tag(n, n);
    bool identity = true;
    for (std::size_t i = 0; i < n; ++i) {
        const int tag = system.atom_tag(i);
        if (tag < 0 || static_cast<std::size_t>(tag) >= n || by_tag[tag] != n) {
            identity = false;
            break;
        }
        by_tag[static_cast<std::size_t>(tag)] = i;
    }
    for (std::size_t i = 0; i < n; ++i) {
        order[i] = identity ? by_tag[i] : i;
    }

    const auto coordinates = system.coordinates();
    const auto forces = system.forces();
    const auto velocities = system.velocities();
    const auto atomic_numbers = system.atomic_numbers();
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = order[k];
        out += std::to_string(system.atom_tag(i));
        out.push_back(' ');
        out += element_symbol(i < atomic_numbers.size() ? atomic_numbers[i] : 0);
        for (int d = 0; d < 3; ++d) {
            out.push_back(' ');
            append_real(out, coordinates[i][d]);
        }
        for (int d = 0; d < 3; ++d) {
            out.push_back(' ');
            append_real(out, i < forces.size() ? forces[i][d] : 0.0);
        }
        if (write_velocities) {
            for (int d = 0; d < 3; ++d) {
                out.push_back(' ');
                // Internal velocity unit is A per internal time unit.
                append_real(out, velocities[i][d] * kInternalTimePerFemtosecond);
            }
        }
        out.push_back('\n');
    }
    return out;
}

void FrameWriter::write_frame(const System& system, std::uint64_t step, double time_fs,
                              double twice_ke, std::size_t dof) {
    if (!out_.is_open()) {
        throw std::runtime_error("FrameWriter: write_frame on a closed writer");
    }
    const std::string frame =
        format_frame(system, step, time_fs, twice_ke, dof, write_velocities_, provenance_);
    out_.write(frame.data(), static_cast<std::streamsize>(frame.size()));
    out_.flush();
    if (!out_) {
        throw std::runtime_error("FrameWriter: failed to write frame at step " +
                                 std::to_string(step));
    }
    ++frame_count_;
}

}  // namespace gmd
