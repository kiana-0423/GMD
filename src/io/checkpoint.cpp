#include "gmd/io/checkpoint.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "gmd/system/system.hpp"

namespace gmd {

namespace {

std::string require_key(std::istream& input, const char* expected) {
    std::string key;
    if (!(input >> key)) {
        throw std::runtime_error(std::string("Checkpoint is missing key '") + expected + "'");
    }
    if (key != expected) {
        throw std::runtime_error(
            "Checkpoint expected key '" + std::string(expected) + "' but found '" + key + "'");
    }
    return key;
}

std::string read_rest_of_line(std::istream& input) {
    std::string value;
    std::getline(input, value);
    if (!value.empty() && value.front() == ' ') {
        value.erase(value.begin());
    }
    return value;
}

void write_topology_section(std::ostream& out, const Topology* topology) {
    if (topology == nullptr) {
        out << "topology_counts 0 0 0 0 0\n";
        return;
    }

    out << "topology_counts "
        << topology->bonds.size() << ' '
        << topology->angles.size() << ' '
        << topology->dihedrals.size() << ' '
        << topology->impropers.size() << ' '
        << topology->constraints.size() << '\n';
    for (const auto& bond : topology->bonds) {
        out << "bond " << bond.i << ' ' << bond.j << ' ' << bond.type_idx << '\n';
    }
    for (const auto& angle : topology->angles) {
        out << "angle " << angle.i << ' ' << angle.j << ' '
            << angle.k << ' ' << angle.type_idx << '\n';
    }
    for (const auto& dihedral : topology->dihedrals) {
        out << "dihedral " << dihedral.i << ' ' << dihedral.j << ' '
            << dihedral.k << ' ' << dihedral.l << ' ' << dihedral.type_idx << '\n';
    }
    for (const auto& improper : topology->impropers) {
        out << "improper " << improper.i << ' ' << improper.j << ' '
            << improper.k << ' ' << improper.l << ' ' << improper.type_idx << '\n';
    }
    for (const auto& constraint : topology->constraints) {
        out << "constraint " << constraint.i << ' ' << constraint.j << ' '
            << constraint.target_distance << '\n';
    }
}

void read_topology_section(std::istream& input, Topology* topology) {
    std::size_t nb = 0;
    std::size_t na = 0;
    std::size_t nd = 0;
    std::size_t ni = 0;
    std::size_t nc = 0;
    require_key(input, "topology_counts");
    if (!(input >> nb >> na >> nd >> ni >> nc)) {
        throw std::runtime_error("Invalid checkpoint topology_counts section");
    }
    if (topology == nullptr) {
        std::string discard;
        std::getline(input, discard);
        for (std::size_t i = 0; i < nb + na + nd + ni + nc; ++i) {
            if (!std::getline(input, discard)) {
                throw std::runtime_error("Checkpoint ended inside topology section");
            }
        }
        return;
    }

    topology->bonds.clear();
    topology->angles.clear();
    topology->dihedrals.clear();
    topology->impropers.clear();
    topology->constraints.clear();

    std::string key;
    for (std::size_t i = 0; i < nb; ++i) {
        BondTerm term;
        if (!(input >> key >> term.i >> term.j >> term.type_idx) || key != "bond") {
            throw std::runtime_error("Invalid checkpoint bond topology entry");
        }
        topology->bonds.push_back(term);
    }
    for (std::size_t i = 0; i < na; ++i) {
        AngleTerm term;
        if (!(input >> key >> term.i >> term.j >> term.k >> term.type_idx) || key != "angle") {
            throw std::runtime_error("Invalid checkpoint angle topology entry");
        }
        topology->angles.push_back(term);
    }
    for (std::size_t i = 0; i < nd; ++i) {
        DihedralTerm term;
        if (!(input >> key >> term.i >> term.j >> term.k >> term.l >> term.type_idx) ||
            key != "dihedral") {
            throw std::runtime_error("Invalid checkpoint dihedral topology entry");
        }
        topology->dihedrals.push_back(term);
    }
    for (std::size_t i = 0; i < ni; ++i) {
        ImproperTerm term;
        if (!(input >> key >> term.i >> term.j >> term.k >> term.l >> term.type_idx) ||
            key != "improper") {
            throw std::runtime_error("Invalid checkpoint improper topology entry");
        }
        topology->impropers.push_back(term);
    }
    for (std::size_t i = 0; i < nc; ++i) {
        BondConstraint constraint;
        if (!(input >> key >> constraint.i >> constraint.j >> constraint.target_distance) ||
            key != "constraint") {
            throw std::runtime_error("Invalid checkpoint constraint topology entry");
        }
        topology->constraints.push_back(constraint);
    }
}

}  // namespace

void write_checkpoint(const std::filesystem::path& path,
                      const CheckpointData& checkpoint) {
    if (checkpoint.system == nullptr) {
        throw std::runtime_error("write_checkpoint requires a System");
    }

    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open checkpoint for writing: " + path.string());
    }
    out << std::setprecision(17);

    const System& system = *checkpoint.system;
    const auto& md = checkpoint.metadata;
    out << "GMD_CHECKPOINT " << kCheckpointVersion << '\n';
    out << "step " << md.step << '\n';
    out << "time_fs " << md.time_fs << '\n';
    out << "boundary " << md.boundary << '\n';
    out << "box " << system.box().lengths[0] << ' '
        << system.box().lengths[1] << ' '
        << system.box().lengths[2] << '\n';
    out << "xyz_file " << md.xyz_file << '\n';
    out << "run_file " << md.run_file << '\n';
    out << "force_field_file " << md.force_field_file << '\n';
    out << "topology_file " << md.topology_file << '\n';
    out << "velocity_seed " << md.velocity_seed << '\n';
    out << "force_field_summary " << md.force_field_summary << '\n';
    out << "config_summary " << md.config_summary << '\n';
    out << "thermostat_type " << md.thermostat_type << '\n';
    out << "thermostat_state " << md.thermostat_state << '\n';
    out << "barostat_type " << md.barostat_type << '\n';
    out << "barostat_state " << md.barostat_state << '\n';
    // Constraint virial state; see CheckpointMetadata for what it means.
    out << "constraint_virial " << md.constraint_virial_state
        << ' ' << md.constraint_virial_time_level;
    for (double component : md.constraint_virial) {
        out << ' ' << component;
    }
    out << '\n';
    out << "atom_count " << system.atom_count() << '\n';
    out << "atoms tag type molecule mass charge x y z vx vy vz\n";

    const auto atom_types = system.atom_types();
    const auto molecules = system.molecule_ids();
    const auto masses = system.masses();
    const auto charges = system.charges();
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();
    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        out << system.atom_tag(atom_index) << ' '
            << atom_types[atom_index] << ' '
            << molecules[atom_index] << ' '
            << masses[atom_index] << ' '
            << charges[atom_index] << ' '
            << coordinates[atom_index][0] << ' '
            << coordinates[atom_index][1] << ' '
            << coordinates[atom_index][2] << ' '
            << velocities[atom_index][0] << ' '
            << velocities[atom_index][1] << ' '
            << velocities[atom_index][2] << '\n';
    }

    write_topology_section(out, checkpoint.topology);
    out << "end\n";
}

CheckpointMetadata read_checkpoint(const std::filesystem::path& path,
                                   System& system,
                                   Topology* topology) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open checkpoint for reading: " + path.string());
    }

    std::string magic;
    int version = 0;
    if (!(input >> magic >> version) || magic != "GMD_CHECKPOINT") {
        throw std::runtime_error("Invalid checkpoint header in " + path.string());
    }
    if (version < kMinReadableCheckpointVersion || version > kCheckpointVersion) {
        throw std::runtime_error(
            "Unsupported checkpoint version " + std::to_string(version) +
            " in " + path.string() + "; this build reads versions " +
            std::to_string(kMinReadableCheckpointVersion) + " to " +
            std::to_string(kCheckpointVersion));
    }

    CheckpointMetadata md;
    require_key(input, "step");
    input >> md.step;
    require_key(input, "time_fs");
    input >> md.time_fs;
    require_key(input, "boundary");
    md.boundary = read_rest_of_line(input);

    require_key(input, "box");
    std::array<double, 3> lengths{};
    if (!(input >> lengths[0] >> lengths[1] >> lengths[2])) {
        throw std::runtime_error("Invalid checkpoint box entry");
    }
    Box box;
    box.set_lengths(lengths);
    require_key(input, "xyz_file");
    md.xyz_file = read_rest_of_line(input);
    require_key(input, "run_file");
    md.run_file = read_rest_of_line(input);
    require_key(input, "force_field_file");
    md.force_field_file = read_rest_of_line(input);
    require_key(input, "topology_file");
    md.topology_file = read_rest_of_line(input);
    require_key(input, "velocity_seed");
    input >> md.velocity_seed;
    require_key(input, "force_field_summary");
    md.force_field_summary = read_rest_of_line(input);
    require_key(input, "config_summary");
    md.config_summary = read_rest_of_line(input);
    require_key(input, "thermostat_type");
    md.thermostat_type = read_rest_of_line(input);
    require_key(input, "thermostat_state");
    md.thermostat_state = read_rest_of_line(input);
    require_key(input, "barostat_type");
    md.barostat_type = read_rest_of_line(input);
    require_key(input, "barostat_state");
    md.barostat_state = read_rest_of_line(input);

    // Version 1 predates the constraint virial state. Such a checkpoint restarts
    // exactly as it always did; its first reported frame simply has no
    // constraint contribution to restore.
    if (version >= 2) {
        require_key(input, "constraint_virial");
        if (!(input >> md.constraint_virial_state >> md.constraint_virial_time_level)) {
            throw std::runtime_error("Invalid checkpoint constraint_virial state");
        }
        if (md.constraint_virial_state != "not_applicable" &&
            md.constraint_virial_state != "unavailable" &&
            md.constraint_virial_state != "valid") {
            throw std::runtime_error("Unknown checkpoint constraint_virial state '" +
                                     md.constraint_virial_state + "'");
        }
        for (double& component : md.constraint_virial) {
            if (!(input >> component)) {
                throw std::runtime_error("Invalid checkpoint constraint_virial entry");
            }
        }

    }

    std::size_t atom_count = 0;
    require_key(input, "atom_count");
    if (!(input >> atom_count)) {
        throw std::runtime_error("Invalid checkpoint atom_count");
    }
    std::string header;
    std::getline(input, header);
    if (!std::getline(input, header) ||
        header != "atoms tag type molecule mass charge x y z vx vy vz") {
        throw std::runtime_error("Invalid checkpoint atoms header");
    }

    system.resize(atom_count);
    system.set_box(box);
    auto atom_types = system.mutable_atom_types();
    auto molecule_ids = system.mutable_molecule_ids();
    auto masses = system.mutable_masses();
    auto charges = system.mutable_charges();
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();
    auto tags = system.mutable_atom_tags();
    auto owners = system.mutable_atom_owners();
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        if (!(input >> tags[atom_index]
              >> atom_types[atom_index]
              >> molecule_ids[atom_index]
              >> masses[atom_index]
              >> charges[atom_index]
              >> coordinates[atom_index][0]
              >> coordinates[atom_index][1]
              >> coordinates[atom_index][2]
              >> velocities[atom_index][0]
              >> velocities[atom_index][1]
              >> velocities[atom_index][2])) {
            throw std::runtime_error("Invalid atom entry in checkpoint");
        }
        owners[atom_index] = 0;
    }

    read_topology_section(input, topology);
    require_key(input, "end");
    system.mutable_neighbor_list().clear();
    return md;
}

}  // namespace gmd
