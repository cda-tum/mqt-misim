#include "dd/MDDPackage.hpp"

#include <chrono>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

dd::Edge<dd::MDDPackage::vNode> fullMixWState(std::ofstream& file, std::vector<size_t> orderOfLayers) {
    std::vector<std::size_t>                      lines{};
    unsigned int                                  numLines = 0U;
    std::map<unsigned int, std::vector<int>>      application;
    std::vector<std::vector<dd::QuantumRegister>> indexes{{}};
    auto                                          sizeTracker = 0;
    auto                                          numop       = 0;
    std::map<unsigned int, unsigned int>          quditcounter;
    quditcounter[2] = 0;
    quditcounter[3] = 0;
    quditcounter[5] = 0;

    bool initial = true;
    for (auto i = 0U; i < orderOfLayers.size(); i++) {
        // Update cardinality of each line, put them in order
        // Update the application indexes
        if (initial) {
            for (auto j = 0U; j < orderOfLayers.at(i); j++) {
                lines.push_back(orderOfLayers.at(i));
                indexes.at(0).push_back(static_cast<dd::QuantumRegister>(j));
                numLines++;
            }
            application[i]     = std::vector<int>{sizeTracker};
            application[i + 1] = {};
            sizeTracker++;
            for (auto j = 0U; j < indexes.at(0).size(); j++) {
                indexes.push_back(std::vector<dd::QuantumRegister>{indexes.at(0).at(static_cast<dd::QuantumRegister>(j))});
                application[i + 1].push_back(sizeTracker);
                sizeTracker++;
            }

            initial = false;
        } else {
            auto tempLine = lines.size();
            auto counter  = 0U;
            for (auto k = 0U; k < tempLine; k++) {
                auto adder = k + counter;
                for (auto j = 1U; j < orderOfLayers.at(i); j++) {
                    lines.insert(lines.begin() + adder + j, orderOfLayers.at(i));
                    numLines++;
                    counter++;
                }
            }
            for (auto c = 0; c < indexes.size(); c++) {
                for (auto f = 0; f < indexes.at(c).size(); f++) {
                    indexes.at(c).at(f) = indexes.at(c).at(f) + indexes.at(c).at(f) * (static_cast<dd::QuantumRegister>(orderOfLayers.at(i) - 1));
                }
            }
            std::vector<dd::QuantumRegister> toAdd{};
            for (auto c = 0; c < indexes.size(); c++) {
                if (indexes.at(c).size() == 1) {
                    toAdd.push_back(indexes.at(c).at(0));
                    for (auto j = 1; j < orderOfLayers.at(i); j++) {
                        indexes.at(c).push_back(indexes.at(c).at(0) + static_cast<dd::QuantumRegister>(j));
                        toAdd.push_back(indexes.at(c).at(0) + static_cast<dd::QuantumRegister>(j));
                    }
                }
            }
            if (i < orderOfLayers.size() - 1) {
                for (auto l = 0; l < toAdd.size(); l++) {
                    indexes.push_back(std::vector<dd::QuantumRegister>{toAdd.at(static_cast<dd::QuantumRegister>(l))});
                    application[i + 1].push_back(sizeTracker);
                    sizeTracker++;
                }
                application[i + 2] = {};
            }
        }
    }

    auto dd = std::make_unique<dd::MDDPackage>(numLines, lines);

    std::vector<size_t> initState(numLines, 0);
    initState.at(0) = 1;

    auto begin = std::chrono::high_resolution_clock::now();

    auto evolution = dd->makeBasisState(numLines, initState);

    for (auto i = 0U; i < orderOfLayers.size(); i++) {
        if (orderOfLayers.at(i) == 2) {
            for (auto g = 0; g < application[i].size(); g++) {
                std::vector<dd::QuantumRegister> inputLines = indexes.at(application[i].at(g));
                evolution                                   = dd->spread2(numLines, inputLines, evolution);
                numop                                       = numop + 3;
            }
        } else if (orderOfLayers.at(i) == 3) {
            for (auto g = 0; g < application[i].size(); g++) {
                evolution = dd->spread3(numLines, reinterpret_cast<const std::vector<dd::QuantumRegister>&>(indexes.at(application[i].at(g))), evolution);
                numop     = numop + 6;
            }
        } else if (orderOfLayers.at(i) == 5) {
            for (auto g = 0; g < application[i].size(); g++) {
                evolution = dd->spread5(numLines, reinterpret_cast<const std::vector<dd::QuantumRegister>&>(indexes.at(application[i].at(g))), evolution);
                numop     = numop + 15;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    std::unordered_set<decltype(evolution.nextNode)> nodeset{};
    auto                                             numnodes = dd->nodeCount(evolution, nodeset);
    auto                                             numcplx  = dd->complexNumber.complexTable.getPeakCount();

    // Build a comma-separated string of the elements
    std::ostringstream oss;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        oss << lines[i];
    }

    for (std::size_t i = 0; i < lines.size(); ++i) {
        quditcounter[lines.at(i)] = quditcounter[lines.at(i)] + 1;
    }

    std::cout << "FullMix, ";
    // Iterate over the map and print the key-value pairs
    for (auto const& pair: quditcounter) {
        std::cout << pair.first << ": " << pair.second << " -";
    }
    std::cout << "//";
    std::cout << lines.size() << ", " << numop << ", " << numnodes << ", " << numcplx << ", " << (elapsed.count() * 1e-9) << "\n";
    return evolution;
}

dd::Edge<dd::MDDPackage::vNode> GHZQutritStateScaled(std::ofstream& file, unsigned int i0) {
    for (auto i = i0; i < i0 + 1; i++) {
        auto                     test = 0;
        std::vector<std::size_t> init(i, 3);
        auto                     dd = std::make_unique<dd::MDDPackage>(i, init);

        auto begin = std::chrono::high_resolution_clock::now();
        // Gates
        auto                               h3Gate = dd->makeGateDD<dd::TritMatrix>(dd::H3(), i, 0);
        std::vector<dd::MDDPackage::mEdge> gates  = {};

        for (auto target = 1U; target < i; target++) {
            dd::Controls target1{};
            dd::Controls target2{};

            for (auto control = 0U; control < target; control++) {
                const dd::Control c1{static_cast<dd::QuantumRegister>(control), 1};
                const dd::Control c2{static_cast<dd::QuantumRegister>(control), 2};
                target1.insert(c1);
                target2.insert(c2);
            }

            gates.push_back(
                    dd->makeGateDD<dd::TritMatrix>(dd::X3, i, target1, target));
            gates.push_back(
                    dd->makeGateDD<dd::TritMatrix>(dd::X3dag, i, target2, target));
        }

        auto evolution = dd->makeZeroState(i);
        evolution      = dd->multiply(h3Gate, evolution);

        for (auto& gate: gates) {
            evolution = dd->multiply(gate, evolution);
        }
        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

        auto                                             numop = gates.size();
        std::unordered_set<decltype(evolution.nextNode)> nodeset{};
        auto                                             numnodes = dd->nodeCount(evolution, nodeset);
        auto                                             numcplx  = dd->complexNumber.complexTable.getPeakCount();

        // Build a comma-separated string of the elements
        std::ostringstream oss;
        for (std::size_t i = 0; i < init.size(); ++i) {
            oss << init[i];
        }

        file << "GHZ, " << init.size() << ", " << oss.str() << ", " << numop << ", " << numnodes << ", " << numcplx << ", " << (elapsed.count() * 1e-9) << "\n";

        return evolution;
    }
}

dd::Edge<dd::MDDPackage::vNode> RandomCircuits(int w, unsigned int d, std::ofstream& file) {
    int          width = w;
    unsigned int depth = d;
    size_t       maxD  = 5;

    std::random_device rd;              // obtain a random number from hardware
    std::mt19937       gen(1592645427); // seed the generator

    std::vector<std::size_t> particles = {};

    std::uniform_int_distribution<> dimdistr(2, maxD);
    for (auto i = 0; i < width; i++) {
        particles.push_back(dimdistr(gen));
    }

    auto dd = std::make_unique<dd::MDDPackage>(width, particles);

    std::uniform_int_distribution<>  pickbool(0, 1);
    std::uniform_int_distribution<>  pickcontrols(1, width - 1);
    std::uniform_real_distribution<> angles(0.0, 2. * dd::PI);

    auto beginClock = std::chrono::high_resolution_clock::now();

    auto evolution = dd->makeZeroState(width);

    for (auto timeStep = 0U; timeStep < depth; timeStep++) {
        for (int line = 0; line < width; line++) {
            //chose if local gate or entangling gate
            auto randomChoice = pickbool(gen);

            if (randomChoice == 0) { //local op

                auto localChoice = pickbool(gen);

                if (localChoice == 0) { //hadamard
                    if (particles.at(line) == 2) {
                        auto chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::H(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::H3(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::H4(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::H5(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                } else { //givens
                    if (particles.at(line) == 2) {
                        long double theta      = 0.;
                        long double phi        = 0.;
                        auto        chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::RXY(theta, phi), width, line);
                        evolution              = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        long double                     theta = angles(gen);
                        long double                     phi   = angles(gen);
                        std::uniform_int_distribution<> picklevel(0, 2);

                        int levelA = picklevel(gen);
                        int levelB = (levelA + 1) % 3;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::RXY3(theta, phi, levelA, levelB), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        long double                     theta = angles(gen);
                        long double                     phi   = angles(gen);
                        std::uniform_int_distribution<> picklevel(0, 3);

                        int levelA = picklevel(gen);
                        int levelB = (levelA + 1) % 4;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::RXY4(theta, phi, levelA, levelB), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        long double                     theta = angles(gen);
                        long double                     phi   = angles(gen);
                        std::uniform_int_distribution<> picklevel(0, 4);

                        int levelA = picklevel(gen);
                        int levelB = (levelA + 1) % 5;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::RXY5(theta, phi, levelA, levelB), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                }
            } else { //entangling

                auto entChoice        = pickbool(gen);
                auto numberOfControls = pickcontrols(gen);

                std::vector<int> controlLines;

                for (auto i = 0U; i < width; i++) {
                    if (i != line) controlLines.push_back(i);
                }

                std::shuffle(begin(controlLines), end(controlLines), gen);
                std::vector<int> controlParticles(controlLines.begin(), controlLines.begin() + numberOfControls);
                std::sort(controlParticles.begin(), controlParticles.end());

                dd::Controls control{};
                for (auto i = 0U; i < numberOfControls; i++) {
                    std::uniform_int_distribution<> picklevel(0, particles.at(controlParticles.at(i)) - 1);
                    auto                            level = picklevel(gen);

                    const dd::Control c{static_cast<dd::QuantumRegister>(controlParticles.at(i)), static_cast<dd::Control::Type>(level)};
                    control.insert(c);
                }

                if (entChoice == 0) { // CEX based
                    //selection of controls
                    if (particles.at(line) == 2) {
                        long double theta      = angles(gen);
                        long double phi        = angles(gen);
                        auto        chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::RXY(theta, phi), width, control, line);
                        evolution              = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        std::uniform_int_distribution<> picklevel(0, 2);

                        long double theta  = angles(gen);
                        long double phi    = angles(gen);
                        int         levelA = picklevel(gen);
                        int         levelB = (levelA + 1) % 3;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::RXY3(theta, phi, levelA, levelB), width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        std::uniform_int_distribution<> picklevel(0, 3);

                        long double theta  = angles(gen);
                        long double phi    = angles(gen);
                        int         levelA = picklevel(gen);
                        int         levelB = (levelA + 1) % 4;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::RXY4(theta, phi, levelA, levelB), width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        std::uniform_int_distribution<> picklevel(0, 4);

                        long double theta  = angles(gen);
                        long double phi    = angles(gen);
                        int         levelA = picklevel(gen);
                        int         levelB = (levelA + 1) % 5;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::RXY5(theta, phi, levelA, levelB), width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                } else { // Controlled clifford
                    if (particles.at(line) == 2) {
                        auto chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::Xmat, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::X3, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::X4, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::X5, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                }
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - beginClock);

    auto                                             numop = width * depth;
    std::unordered_set<decltype(evolution.nextNode)> nodeset{};
    auto                                             numnodes = dd->nodeCount(evolution, nodeset);
    auto                                             numcplx  = dd->complexNumber.complexTable.getPeakCount();

    // Build a comma-separated string of the elements
    std::ostringstream oss;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        oss << particles[i];
    }

    file << "Random, " << particles.size() << ", " << oss.str() << ", " << numop << ", " << numnodes << ", " << numcplx << ", " << (elapsed.count() * 1e-9) << "\n";
    return evolution;
}
int main() {
    std::ofstream myfile;
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);
    myfile << "Bench, NumLines, Qudits, Operations, Nodes, CplxNum, time\n";

    fullMixWState(myfile, {2, 3});
    fullMixWState(myfile, {3, 5, 2});
    fullMixWState(myfile, {3, 5, 2, 3});
    fullMixWState(myfile, {2, 3, 2, 3, 3});
    fullMixWState(myfile, {2, 5, 3, 3});
    fullMixWState(myfile, {2, 3, 2, 5});
    fullMixWState(myfile, {2, 3, 3, 3});
    fullMixWState(myfile, {5, 5, 2, 2});
    std::cout << "First set"
              << "\n";
    myfile.close();
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);

    GHZQutritStateScaled(myfile, 5);
    GHZQutritStateScaled(myfile, 10);
    GHZQutritStateScaled(myfile, 30);
    GHZQutritStateScaled(myfile, 60);
    GHZQutritStateScaled(myfile, 120);
    GHZQutritStateScaled(myfile, 128);
    std::cout << "Second set"
              << "\n";
    myfile.close();
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);
    RandomCircuits(3, 1000, myfile);
    RandomCircuits(4, 1000, myfile);
    RandomCircuits(5, 1000, myfile);
    RandomCircuits(6, 1000, myfile);
    RandomCircuits(7, 1000, myfile);

    std::cout << "Third set"
              << "\n";
    myfile.close();
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);
    RandomCircuits(3, 1000, myfile);
    RandomCircuits(4, 1000, myfile);
    RandomCircuits(5, 1000, myfile);
    RandomCircuits(6, 1000, myfile);
    RandomCircuits(7, 1000, myfile);

    std::cout << "fourth set"
              << "\n";
    myfile.close();
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);
    RandomCircuits(3, 1000, myfile);
    RandomCircuits(4, 1000, myfile);
    RandomCircuits(5, 1000, myfile);
    RandomCircuits(6, 1000, myfile);
    RandomCircuits(7, 1000, myfile);

    std::cout << "fifth set"
              << "\n";
    myfile.close();
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);
    RandomCircuits(8, 1000, myfile);
    RandomCircuits(9, 1000, myfile);
    RandomCircuits(10, 1000, myfile);
    RandomCircuits(11, 1000, myfile);
    RandomCircuits(12, 1000, myfile);
    RandomCircuits(13, 1000, myfile);

    std::cout << "6 set"
              << "\n";
    myfile.close();
    myfile.open("/home/k3vn/Desktop/trycollect.csv", std::ios_base::app);
    RandomCircuits(8, 1000, myfile);
    RandomCircuits(9, 1000, myfile);
    RandomCircuits(10, 1000, myfile);
    RandomCircuits(11, 1000, myfile);
    RandomCircuits(12, 1000, myfile);
    RandomCircuits(13, 1000, myfile);

    std::cout << "7 set"
              << "\n";

    myfile.close();
    return 0;
}
