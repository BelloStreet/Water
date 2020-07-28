#include "../include/MoldenParser.hpp"
#include <boost/algorithm/string_regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <cstring>
#include <fstream>
#include <iostream>

/**
  Function that parses the Molden input file and populates the global datatypes
  molecule_t and MO_t. This was determined to be a function because it "does"
  something (i.e. a function) and a class "is" something.
  */
void moldenParser(const std::string &a_filename, molecule_t &a_molecule) {

  /* TODO:should parse the units */
  /* double angstrom_flag = 1.0; */
  double angstrom_flag = 1.889726124564062;
  /* if (a_filename.find("Angs") != std::string::npos) { */
  /*   std::cout << " !!!found angs " */
  /*             << "\n"; */
  /* } */

  // Seperating file into strings
  boost::filesystem::ifstream source(a_filename);
  std::string MolproString, ParticlesString, GaussianString, MOString;
  std::string line;
  while (std::getline(source, line)) {
    if (line.find("Molpro") != std::string::npos) {
      std::getline(source, line, '['); // reads until [
      MolproString.append(line);
    } else if (line.find("Atoms") != std::string::npos) {
      std::getline(source, line, '['); // reads until [
      ParticlesString.append(line);
    } else if (line.find("GTO") != std::string::npos) {
      std::getline(source, line, '['); // reads until [
      GaussianString.append(line);
    } else if (line.find("MO") != std::string::npos) {
      std::getline(source, line, '\0'); // reads until end of file
      MOString.append(line);
    }
  }

  /* Getting Point Group if possible */
  if (MolproString.find("Cs") != std::string::npos) {
    a_molecule.pointgroup = "Cs";
  } else if (MolproString.find("C2v") != std::string::npos) {
    a_molecule.pointgroup = "C2v";
  } else if (MolproString.find("D2h") != std::string::npos) {
    a_molecule.pointgroup = "D2h";
  } else {
    std::cout
        << " Couldn't find point group!  Need to see if the "
           "molecular orbitals' symmetries have the conventional symbols.  "
           "If "
           "not, then may need to look up convention used by the quantum "
           "packaged that generated the Molden file, like Molpro (see "
           "https://www.molpro.net/info/2015/doc/manual/node36.html) "
        << "\n";
  }

  /* Saving atom information */
  std::vector<std::string> symbol;
  std::vector<int> charge;
  std::vector<double> x, y, z;
  boost::char_separator<char> new_line{"\n"};
  boost::tokenizer<boost::char_separator<char>> particle_tok{ParticlesString,
                                                             new_line};
  for (auto full_line = particle_tok.begin(); full_line != particle_tok.end();
       ++full_line) {
    std::vector<std::string> atoms;
    boost::split(atoms, *full_line, boost::is_any_of(" "));
    atoms.erase(std::remove(atoms.begin(), atoms.end(), ""), atoms.end());
    symbol.push_back(atoms[0]);
    charge.push_back(std::stoi(atoms[2]));
    x.push_back(std::stod(atoms[3]));
    y.push_back(std::stod(atoms[4]));
    z.push_back(std::stod(atoms[5]));
  }

  int number_atoms = symbol.size();
  molecule_t molecule;
  molecule.atoms.reserve(number_atoms);

  /* Saving Gaussian information */
  int nbasis;
  std::vector<std::string> particles;
  std::vector<std::string> token;
  auto re_blank_space = boost::regex(R"(\n\s*\n)");
  boost::split_regex(particles, GaussianString, re_blank_space);
  atom_t tmp_atom;
  for (int particle = 0; particle < number_atoms; ++particle) {
    boost::split(token, particles[particle], boost::is_any_of(" "));
    token.erase(std::remove(token.begin(), token.end(), ""), token.end());
    AOorbital_t tmp_AOorbital;
    for (int i = 0; i < token.size(); ++i) {
      if (token[i] == "s") {
        Sorbitalparameters_t gaussians;
        nbasis = nbasis + 1;
        std::vector<double> tmp;
        tmp.reserve(std::stoi(token[i + 1]));
        gaussians.expon_s.reserve(std::stoi(token[i + 1]));
        gaussians.coeff_s.reserve(std::stoi(token[i + 1]));
        /* Strange loop mapping is to stride the tokens correctly */
        for (int j = i + 3; j < i + (2 * std::stoi(token[i + 1]) + 3); ++j) {
          std::replace(token[j].begin(), token[j].end(), 'D', 'E');
          tmp.push_back(std::stod(token[j]));
        }
        for (int j = 0; j < tmp.size(); ++j) {
          if (j % 2 == 0) {
            gaussians.expon_s.push_back(tmp[j]);
          } else {
            auto fact = pow(2.0 * tmp[j - 1] / M_PI, 0.75);
            gaussians.coeff_s.push_back(fact * tmp[j]);
          }
        }
        tmp_AOorbital.s_orbitals.push_back(gaussians);
      }
      if (token[i] == "p") {
        Porbitalparameters_t gaussians;
        nbasis = nbasis + 3;
        std::vector<double> tmp;
        tmp.reserve(std::stoi(token[i + 1]));
        gaussians.expon_p.reserve(std::stoi(token[i + 1]));
        gaussians.coeff_p.reserve(std::stoi(token[i + 1]));
        for (int j = i + 3; j < i + (2 * std::stoi(token[i + 1]) + 3); ++j) {
          std::replace(token[j].begin(), token[j].end(), 'D', 'E');
          tmp.push_back(std::stod(token[j]));
        }
        for (int j = 0; j < tmp.size(); ++j) {
          if (j % 2 == 0) {
            gaussians.expon_p.push_back(tmp[j]);
          } else {
            auto fact = pow(
                pow(2.0, 7.0) * pow(tmp[j - 1], 5.0) / pow(M_PI, 3.0), 0.25);
            gaussians.coeff_p.push_back(fact * tmp[j]);
          }
        }
        tmp_AOorbital.p_orbitals.push_back(gaussians);
      }
      if (token[i] == "d") {
        Dorbitalparameters_t gaussians;
        nbasis = nbasis + 6;
        std::vector<double> tmp;
        tmp.reserve(std::stoi(token[i + 1]));
        for (int j = i + 3; j < i + (2 * std::stoi(token[i + 1]) + 3); ++j) {
          std::replace(token[j].begin(), token[j].end(), 'D', 'E');
          tmp.push_back(std::stod(token[j]));
        }
        for (int j = 0; j < tmp.size(); ++j) {
          if (j % 2 == 0) {
            gaussians.expon_d.push_back(tmp[j]);
          } else {
            auto fact = pow(
                pow(2.0, 11.0) * pow(tmp[j - 1], 7.0) / pow(M_PI, 3.0), 0.25);
            gaussians.coeff_d.push_back(fact * tmp[j]);
          }
        }
        tmp_AOorbital.d_orbitals.push_back(gaussians);
      }
      if (token[i] == "f") {
        Forbitalparameters_t gaussians;
        nbasis = nbasis + 10;
        std::vector<double> tmp;
        tmp.reserve(std::stoi(token[i + 1]));
        gaussians.expon_f.reserve(std::stoi(token[i + 1]));
        gaussians.coeff_f.reserve(std::stoi(token[i + 1]));
        for (int j = i + 3; j < i + (2 * std::stoi(token[i + 1]) + 3); ++j) {
          std::replace(token[j].begin(), token[j].end(), 'D', 'E');
          tmp.push_back(std::stod(token[j]));
        }
        for (int j = 0; j < tmp.size(); ++j) {
          if (j % 2 == 0) {
            gaussians.expon_f.push_back(tmp[j]);
          } else {
            auto fact = pow(
                pow(2.0, 15.0) * pow(tmp[j - 1], 9.0) / pow(M_PI, 3.0), 0.25);
            gaussians.coeff_f.push_back(fact * tmp[j]);
          }
        }
        tmp_AOorbital.f_orbitals.push_back(gaussians);
      }
    }
    tmp_atom.AOorbital = tmp_AOorbital;
    tmp_atom.symbol = symbol[particle];
    tmp_atom.x = x[particle] * angstrom_flag;
    tmp_atom.y = y[particle] * angstrom_flag;
    tmp_atom.z = z[particle] * angstrom_flag;
    tmp_atom.charge = charge[particle];
    a_molecule.atoms.push_back(tmp_atom);
  }

  /* Saving molecular orbital information */
  std::vector<MO_t> tmp_MOs;
  std::vector<std::string> molecular_orbitals;
  auto re_sym = boost::regex(R"(Sym=)");
  boost::split_regex(molecular_orbitals, MOString, re_sym);
  for (int mo = 1; mo < molecular_orbitals.size(); ++mo) {
    MO_t tmp_mo;
    tmp_mo.value.reserve(nbasis - 1);
    std::vector<std::string> token;
    boost::split(token, molecular_orbitals[mo], boost::is_any_of("\n"));
    token.erase(std::remove(token.begin(), token.end(), ""), token.end());
    std::vector<std::string> occupation;
    boost::split(occupation, token[3], boost::is_any_of(" "));
    occupation.erase(std::remove(occupation.begin(), occupation.end(), ""),
                     occupation.end());
    if (std::atof(occupation[1].c_str()) != 0.0) {
      tmp_mo.occupation = std::atof(occupation[1].c_str());
      std::vector<std::string> symmetry;
      boost::split(symmetry, token[0], boost::is_any_of(" "));
      symmetry.erase(std::remove(symmetry.begin(), symmetry.end(), ""),
                     symmetry.end());
      if (a_molecule.pointgroup == "C2v") {
        if (symmetry[0].find(".1") != std::string::npos) {
          tmp_mo.symmetry = "A1";
        } else if (symmetry[0].find(".2") != std::string::npos) {
          tmp_mo.symmetry = "B1";
        } else if (symmetry[0].find(".3") != std::string::npos) {
          tmp_mo.symmetry = "B2";
        } else if (symmetry[0].find(".4") != std::string::npos) {
          tmp_mo.symmetry = "A2";
        } else {
          std::cout << " Didn't change Molpro convention in C2v symmetry. "
                    << "\n";
        }
      } else if (a_molecule.pointgroup == "D2h") {
        if (symmetry[0].find(".1") != std::string::npos) {
          tmp_mo.symmetry = "Ag";
        } else if (symmetry[0].find(".2") != std::string::npos) {
          tmp_mo.symmetry = "B3u";
        } else if (symmetry[0].find(".3") != std::string::npos) {
          tmp_mo.symmetry = "B2u";
        } else if (symmetry[0].find(".4") != std::string::npos) {
          tmp_mo.symmetry = "B1g";
        } else if (symmetry[0].find(".5") != std::string::npos) {
          tmp_mo.symmetry = "B1u";
        } else if (symmetry[0].find(".6") != std::string::npos) {
          tmp_mo.symmetry = "B2g";
        } else if (symmetry[0].find(".7") != std::string::npos) {
          tmp_mo.symmetry = "B3g";
        } else if (symmetry[0].find(".8") != std::string::npos) {
          tmp_mo.symmetry = "Au";
        } else {
          std::cout << " Didn't change Molpro convention in D2h symmetry. "
                    << "\n";
        }
      } else if (a_molecule.pointgroup == "Cs") {
        if (symmetry[0].find(".1") != std::string::npos) {
          tmp_mo.symmetry = "A'";
        } else if (symmetry[0].find(".2") != std::string::npos) {
          tmp_mo.symmetry = "A''";
        } else {
          std::cout << " Didn't change Molpro convention in Cs symmetry. "
                    << "\n";
        }
      } else {
        std::cout << " Not using Molpro to generate Molden file.  Molecule "
                     "point group will be NAN. "
                  << "\n";
      }
      std::vector<std::string> energy;
      boost::split(energy, token[1], boost::is_any_of(" "));
      energy.erase(std::remove(energy.begin(), energy.end(), ""), energy.end());
      tmp_mo.energy = std::atof(energy[1].c_str());
      std::vector<std::string> spin;
      boost::split(spin, token[2], boost::is_any_of(" "));
      spin.erase(std::remove(spin.begin(), spin.end(), ""), spin.end());
      tmp_mo.spin = spin[1];
    }
    std::vector<std::string> value_token;
    if (std::atof(occupation[1].c_str()) != 0.0) {
      for (int i = 4; i < token.size() - 1; ++i) {
        boost::split(value_token, token[i], boost::is_any_of(" "));
        value_token.erase(
            std::remove(value_token.begin(), value_token.end(), ""),
            value_token.end());
        tmp_mo.value.push_back(std::atof(value_token[1].c_str()));
      }
      tmp_MOs.push_back(tmp_mo);
    }
  }
  a_molecule.molecular_orbitals = tmp_MOs;
}
