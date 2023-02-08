/*
 * This file is part of the MQT DD MDDPackage which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */
#ifndef DDMDDPackage_H
#define DDMDDPackage_H

#include "Complex.hpp"
#include "ComplexTable.hpp"
#include "ComplexValue.hpp"
#include "ComputeTable.hpp"
#include "ComplexNumbers.hpp"
#include "Control.hpp"
#include "Definitions.hpp"
#include "Edge.hpp"
#include "GateMatrixDefinitions.hpp"
#include "UniqueTable.hpp"
#include "UnaryComputeTable.hpp"


#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>


namespace dd {
    class MDDPackage {
  ///
  /// Complex number handling
  ///
 public:
  ComplexNumbers complexNumber{};

  ///
  /// Construction, destruction, information and reset
  ///
 public:
  static constexpr std::size_t MAX_POSSIBLE_REGISTERS =
      static_cast<std::make_unsigned_t<QuantumRegister>>(
          std::numeric_limits<QuantumRegister>::max()) +
      1U;
  static constexpr std::size_t DEFAULT_REGISTERS = 128;

  explicit MDDPackage(std::size_t nqr, std::vector<size_t> sizes)
      : complexNumber(ComplexNumbers()),
        numberOfQuantumRegisters(nqr),
        registersSizes(std::move(sizes)) {
    resize(nqr);
  };

  ~MDDPackage() = default;

  MDDPackage(const MDDPackage& MDDPackage) = delete;  // no copy constructor
  MDDPackage& operator=(const MDDPackage& MDDPackage) =
      delete;  // no copy assignment constructor

  // TODO RESIZE
  //  resize the package instance
  void resize(std::size_t nq) {
    // TODO DISCUSS THIS FEATURE
    if (nq > MAX_POSSIBLE_REGISTERS) {
      throw std::invalid_argument(
          "Requested too many qubits from package. Qubit datatype only "
          "allows up to " +
          std::to_string(MAX_POSSIBLE_REGISTERS) + " qubits, while " +
          std::to_string(nq) +
          " were requested. Please recompile the package with a wider "
          "Qubit type!");
    }
    numberOfQuantumRegisters = nq;
    vUniqueTable.resize(numberOfQuantumRegisters);
    mUniqueTable.resize(numberOfQuantumRegisters);
    // dUniqueTable.resize(number_of_quantum_registers);
    // stochasticNoiseOperationCache.resize(number_of_quantum_registers);
    idTable.resize(numberOfQuantumRegisters);
  }

  // reset package state
  void reset() {
// TODO IMPLEMENT
// clearUniqueTables();
// clearComputeTables();
complexNumber.clear();
  }

  // TODO CHECK SYNTAX OF SETTERS AND GETTERS
  //  getter for number qudits

  [[nodiscard]] auto qregisters() const { return numberOfQuantumRegisters; }

  // setter dimensionalisties
  [[nodiscard]] auto registerDimensions(const std::vector<size_t>& regs) {
    registersSizes = regs;
  }

  // getter for sizes
  [[nodiscard]] auto regsSize() const { return registersSizes; }

 private:
  std::size_t numberOfQuantumRegisters;
  // TODO THIS IS NOT CONST RIGHT?
  // from LSB TO MSB
  std::vector<size_t> registersSizes;

  ///
  /// Vector nodes, edges and quantum states
  ///
 public:
  struct vNode {
    std::vector<Edge<vNode>> edges{};  // edges out of this node
    vNode* next{};                     // used to link nodes in unique table
    RefCount refCount{};  // reference count, how many active dd are using
                          // the node
    QuantumRegister
        varIndx{};  // variable index (nonterminal) value (-1 for terminal),
                    // index in the circuit endianess 0 from below

    static vNode terminalNode;
    constexpr static vNode* terminal{&terminalNode};

    static constexpr bool isTerminal(const vNode* nodePoint) {
      return nodePoint == terminal;
    }
  };
  using vEdge = Edge<vNode>;
  using vCachedEdge = CachedEdge<vNode>;

  vEdge normalize(const vEdge& edge, bool cached) {
    std::vector<bool> zero;

    for (auto& i : edge.nextNode->edges) {
      zero.push_back(i.weight.approximatelyZero());
    }

    // make sure to release cached numbers approximately zero, but not exactly
    // zero
    if (cached) {
      for (auto i = 0U; i < zero.size(); i++) {
        if (zero.at(i) && edge.nextNode->edges.at(i).weight != Complex::zero) {
          // TODO what is returnToCache

          complexNumber.returnToCache(edge.nextNode->edges.at(i).weight);
          edge.nextNode->edges.at(i) = vEdge::zero;
        }
      }
    }

            // all equal to zero
            if ( none_of( cbegin(zero), cend(zero), std::logical_not<bool>() ) ) {

                if (!cached && !edge.isTerminal()) {
                    // If it is not a cached computation, the node has to be put back into the chain
                    vUniqueTable.returnNode(edge.nextNode);
                }
                return vEdge::zero;
            }

            // find indices that are not zero
            std::vector<int> nonZeroIndices;
            for (auto i = 0U; i < zero.size(); i++) {
              if (!zero.at(i)) {
                nonZeroIndices.push_back(i);
              }
            }

            if (nonZeroIndices.size() == 1) {
              // search for first element different from zero
              auto currentEdge = edge;
              auto& weightFromChild =
                  currentEdge.nextNode->edges[nonZeroIndices.front()].weight;

              if (cached && weightFromChild != Complex::one) {
                currentEdge.weight = weightFromChild;
              } else {
                currentEdge.weight = complexNumber.lookup(weightFromChild);
              }

              weightFromChild = Complex::one;
              return currentEdge;
            }

            // calculate normalizing factor
            auto sumNorm2 =
                ComplexNumbers::mag2(edge.nextNode->edges.at(0).weight);

            for (auto i = 1; i < edge.nextNode->edges.size(); i++) {
              sumNorm2 = sumNorm2 + ComplexNumbers::mag2(
                                        edge.nextNode->edges.at(i).weight);
            }

            const auto commonFactor = std::sqrt(sumNorm2);

            // set incoming edge weight to max
            auto currentEdge = edge;

            // TODO CHECK EXACTLY THIS CACHING SYSTEM
            // auto& max = r.nextNode->edges[argMax];
            if (cached && currentEdge.weight != Complex::one) {
              // current_edge.weight = max.weight;
              currentEdge.weight.real->value *= commonFactor;
              currentEdge.weight.img->value *= commonFactor;
            } else {
              currentEdge.weight = complexNumber.lookup(
                  CTEntry::val(currentEdge.weight.real) * commonFactor,
                  CTEntry::val(currentEdge.weight.img) * commonFactor);
              if (currentEdge.weight.approximatelyZero()) {
                return vEdge::zero;
              }
            }

            // actual normalization of the edges
            for (auto& i : edge.nextNode->edges) {
              if (cached) {
                complexNumber.returnToCache(i.weight);
                ComplexNumbers::div(i.weight, i.weight, currentEdge.weight);
                i.weight = complexNumber.lookup(i.weight);
              } else {
                auto c = complexNumber.getTemporary();
                ComplexNumbers::div(c, i.weight, currentEdge.weight);
                i.weight = complexNumber.lookup(c);
              }
              if (i.weight == Complex::zero) {
                i = vEdge::zero;
              }
            }

            return currentEdge;
  }

  // generate |0...0> with N quantum registers
        vEdge makeZeroState(QuantumRegisterCount n, std::size_t start = 0) {
    if (n + start > numberOfQuantumRegisters) {
      // TODO UNDERSTAND RESIZING
      throw std::runtime_error("Requested state with " +
                               std::to_string(n + start) +
                               " QUANTUM REGISTERS, but current package "
                               "configuration only supports up to " +
                               std::to_string(numberOfQuantumRegisters) +
                               " QUANTUM REGISTERS. Please allocate a "
                               "larger package instance.");
    }
    auto first = vEdge::one;
    for (std::size_t node_idx = start; node_idx < n + start; node_idx++) {
      std::vector<Edge<vNode>> new_outgoing_edges;
      new_outgoing_edges.reserve(registersSizes.at(node_idx));
      new_outgoing_edges.push_back(first);
      for (auto i = 1U; i < registersSizes.at(node_idx); i++) {
        new_outgoing_edges.push_back(vEdge::zero);
      }

      first = makeDDNode(static_cast<QuantumRegister>(node_idx),
                         new_outgoing_edges);
    }
    return first;
  }
  // generate computational basis state |i> with n quantum registers
  vEdge makeBasisState(QuantumRegisterCount n, const std::vector<bool>& state,
                       std::size_t start = 0) {
    if (n + start > numberOfQuantumRegisters) {
      throw std::runtime_error(
          "Requested state with " + std::to_string(n + start) +
          " qubits, but current package configuration only supports up "
          "to " +
          std::to_string(numberOfQuantumRegisters) +
          " qubits. Please allocate a larger package instance.");
    }
    auto f = vEdge::one;
    for (std::size_t pos = start; pos < n + start; ++pos) {
      if (state.at(pos) == 0) {
        f = makeDDNode(static_cast<QuantumRegister>(pos),
                       std::vector{f, vEdge::zero});
      } else {
        f = makeDDNode(static_cast<QuantumRegister>(pos),
                       std::vector{vEdge::zero, f});
      }
    }
    return f;
  }

  // create a normalized DD node and return an edge pointing to it. The
  // node is not recreated if it already exists.
  template <class Node>
  Edge<Node> makeDDNode(QuantumRegister varidx,
                        const std::vector<Edge<Node>>& edges,
                        bool cached = false) {
    auto& uniqueTable = getUniqueTable<Node>();

    Edge<Node> new_edge{uniqueTable.getNode(), Complex::one};
    new_edge.nextNode->varIndx = varidx;
    new_edge.nextNode->edges = edges;

    assert(new_edge.nextNode->refCount == 0);

    for ([[maybe_unused]] const auto& edge : edges)
      assert(edge.nextNode->varIndx == varidx - 1 || edge.isTerminal());

    // normalize it
    new_edge = normalize(new_edge, cached);
    assert(new_edge.nextNode->varIndx == varidx || new_edge.isTerminal());

    // look it up in the unique tables
    auto looked_up_edge = uniqueTable.lookup(new_edge, false);
    assert(looked_up_edge.nextNode->varIndx == varidx ||
           looked_up_edge.isTerminal());

    // set specific node properties for matrices
    if constexpr (std::is_same_v<Node, mNode>) {
      if (looked_up_edge.nextNode == new_edge.nextNode) {
      }

      checkSpecialMatrices(looked_up_edge.nextNode);
    }

    return looked_up_edge;
  }

 public:
  struct mNode {
    std::vector<Edge<mNode>> edges{};  // edges out of this node
    mNode* next{};                     // used to link nodes in unique table
    RefCount refCount{};               // reference count
    QuantumRegister varIndx{};         // variable index (nonterminal) value (-1
                                       // for terminal)
    bool symmetric = false;            // node is symmetric
    bool identity = false;             // node resembles identity

    static mNode terminalNode;
    constexpr static mNode* terminal{&terminalNode};

    static constexpr bool isTerminal(const mNode* node_point) {
      return node_point == terminal;
    }
  };
        using mEdge       = Edge<mNode>;
        using mCachedEdge = CachedEdge<mNode>;


        mEdge normalize(const mEdge& edge, bool cached) {
          auto argmax = -1;

          std::vector<bool> zero;

          for (auto& i : edge.nextNode->edges) {
            zero.push_back(i.weight.approximatelyZero());
          }

          // make sure to release cached numbers approximately zero, but not
          // exactly zero
          if (cached) {
            for (auto i = 0U; i < zero.size(); i++) {
              if (zero.at(i) &&
                  edge.nextNode->edges.at(i).weight != Complex::zero) {
                // TODO what is returnToCache

                complexNumber.returnToCache(edge.nextNode->edges.at(i).weight);
                edge.nextNode->edges.at(i) = mEdge::zero;
              }
            }
          }

            fp   max_magnitude  = 0;
            auto max_weight = Complex::one;
            // determine max amplitude
            for (auto i = 0U; i < zero.size(); ++i) {
                if (zero.at(i)) continue;
                if (argmax == -1) {
                  argmax = static_cast<decltype(argmax)>(i);
                  max_magnitude =
                      ComplexNumbers::mag2(edge.nextNode->edges.at(i).weight);
                  max_weight = edge.nextNode->edges.at(i).weight;
                } else {
                  auto current_magnitude =
                      ComplexNumbers::mag2(edge.nextNode->edges.at(i).weight);
                  if (current_magnitude - max_magnitude >
                      ComplexTable<>::tolerance()) {
                    argmax = static_cast<decltype(argmax)>(i);
                    max_magnitude = current_magnitude;
                    max_weight = edge.nextNode->edges.at(i).weight;
                  }
                }
            }

            // all equal to zero
            if (argmax == -1) {
              if (!cached && !edge.isTerminal()) {
                // If it is not a cached computation, the node has to be put
                // back into the chain
                mUniqueTable.returnNode(edge.nextNode);
              }
              return mEdge::zero;
            }

            auto current_edge = edge;
            // divide each entry by max
            for (auto i = 0U; i < edge.nextNode->edges.size(); ++i) {
              if (static_cast<decltype(argmax)>(i) == argmax) {
                if (cached) {
                  if (current_edge.weight == Complex::one)
                    current_edge.weight = max_weight;
                  else
                    ComplexNumbers::mul(current_edge.weight,
                                        current_edge.weight, max_weight);
                } else {
                  if (current_edge.weight == Complex::one) {
                    current_edge.weight = max_weight;
                  } else {
                    auto new_complex_numb = complexNumber.getTemporary();
                    ComplexNumbers::mul(new_complex_numb, current_edge.weight,
                                        max_weight);
                    current_edge.weight =
                        complexNumber.lookup(new_complex_numb);
                  }
                }
                current_edge.nextNode->edges.at(i).weight = Complex::one;
              } else {
                if (cached && !zero.at(i) &&
                    current_edge.nextNode->edges.at(i).weight != Complex::one) {
                  complexNumber.returnToCache(
                      current_edge.nextNode->edges.at(i).weight);
                }
                if (current_edge.nextNode->edges.at(i)
                        .weight.approximatelyOne())
                  current_edge.nextNode->edges.at(i).weight = Complex::one;
                auto new_complex_numb = complexNumber.getTemporary();

                ComplexNumbers::div(new_complex_numb,
                                    current_edge.nextNode->edges.at(i).weight,
                                    max_weight);
                current_edge.nextNode->edges.at(i).weight =
                    complexNumber.lookup(new_complex_numb);
              }
            }
            return current_edge;
        }

        /// Make GATE DD
        // SIZE => EDGE (number of successors)
        // build matrix representation for a single gate on an n-qubit circuit

        template <typename Matrix>
        mEdge makeGateDD(const Matrix& mat, QuantumRegisterCount n,
                         QuantumRegister target, std::size_t start = 0) {
          return makeGateDD(mat, n, Controls{}, target, start);
        }

        template <typename Matrix>
        mEdge makeGateDD(const Matrix& mat, QuantumRegisterCount n,
                         const Control& control, QuantumRegister target,
                         std::size_t start = 0) {
          return makeGateDD(mat, n, Controls{control}, target, start);
        }

        template <typename Matrix>
        mEdge makeGateDD(const Matrix& mat, QuantumRegisterCount n,
                         const Controls& controls, QuantumRegister target,
                         std::size_t start = 0) {
          if (n + start > numberOfQuantumRegisters) {
            throw std::runtime_error(
                "Requested gate with " + std::to_string(n + start) +
                " qubits, but current package configuration only supports up "
                "to " +
                std::to_string(numberOfQuantumRegisters) +
                " qubits. Please allocate a larger package instance.");
          }

          auto targetRadix =
              registersSizes.at(static_cast<unsigned long>(target));
          auto edges = targetRadix * targetRadix;
          std::vector<mEdge> edgesMat{};
          auto currentControl = controls.begin();

          for (auto i = 0U; i < edges; ++i) {
            if (mat.at(i).r == 0 && mat.at(i).i == 0) {
              edgesMat.at(i) = mEdge::zero;
            } else {
              edgesMat.at(i) = mEdge::terminal(complexNumber.lookup(mat.at(i)));
            }
          }
          auto currentReg = static_cast<QuantumRegister>(start);
          // process lines below target
          for (; currentReg < target; currentReg++) {
            auto radix =
                registersSizes.at(static_cast<unsigned long>(currentReg));

            for (auto rowMat = 0U; rowMat < targetRadix; ++rowMat) {
              for (auto colMat = 0U; colMat < targetRadix; ++colMat) {
                auto entryPos = rowMat + targetRadix * colMat;

                std::vector<mEdge> quadEdges(radix * radix, mEdge::zero);

                if (currentControl != controls.end() &&
                    currentControl->quantumRegister == currentReg) {
                  // TODO TAKE CARE OF PATTERNS, OR LIMIT WHERE THEY CAN PUT
                  // CONTROLS

                  if (rowMat == colMat) {
                    for (auto i = 0U; i < radix; i++) {
                      auto diagInd = i * radix + i;
                      if (diagInd == currentControl->type) {
                        quadEdges.at(diagInd) = edgesMat.at(entryPos);
                      } else {
                        quadEdges.at(diagInd) = mEdge::one;
                      }
                    }
                  } else {
                    quadEdges.at(currentControl->type +
                                 radix * currentControl->type) =
                        edgesMat.at(entryPos);
                  }
                  edgesMat.at(entryPos) = makeDDNode(currentReg, quadEdges);

                } else {  // not connected

                  for (auto iD = 0U; iD < radix; iD++) {
                    quadEdges.at(iD * radix + iD) = edgesMat.at(entryPos);
                  }
                  edgesMat.at(entryPos) = makeDDNode(currentReg, quadEdges);
                }
              }
            }

            if (currentControl != controls.end() &&
                currentControl->quantumRegister == currentReg) {
              ++currentReg;
            }
          }

          // target line
          auto targetNodeEdge = makeDDNode(currentReg, edgesMat);

          // process lines above target
          for (; currentReg < static_cast<QuantumRegister>(n - 1 + start);
               currentReg++) {
            auto nextReg = static_cast<QuantumRegister>(currentReg + 1);
            auto nextRadix =
                registersSizes.at(static_cast<unsigned long>(currentReg));
            std::vector<mEdge> nextEdges(nextRadix * nextRadix, mEdge::zero);

            if (currentControl != controls.end() &&
                currentControl->quantumRegister == nextReg) {
              for (auto i = 0U; i < nextRadix; i++) {
                auto diagInd = i * nextRadix + i;
                if (diagInd == currentControl->type) {
                  nextEdges.at(diagInd) = targetNodeEdge;
                } else {
                  nextEdges.at(diagInd) = mEdge::one;
                }
              }
              /*
              if (currentControl->type ) { // neg. control
                targetNodeEdge = makeDDNode(nextReg, std::array{
                        targetNodeEdge, mEdge::zero, mEdge::zero,
              makeIdent(static_cast<QuantumRegister>(start),
              static_cast<QuantumRegister>(nextReg - 1))}); } else { // pos.
              control targetNodeEdge = makeDDNode(nextReg,
              std::array{makeIdent(static_cast<QuantumRegister>(start),
              static_cast<QuantumRegister>(nextReg - 1)), mEdge::zero,
              mEdge::zero, targetNodeEdge});
              }
              */

              ++currentControl;

            } else {  // not connected
              for (auto iD = 0U; iD < nextRadix; iD++) {
                nextEdges.at(iD * nextRadix + iD) = targetNodeEdge;
              }
              targetNodeEdge = makeDDNode(nextReg, nextEdges);
            }
          }
          return targetNodeEdge;
        }

        ///
        /// Identity matrices
        ///
       public:
        // create n-qudit identity DD. makeIdent(n) === makeIdent(0, n-1)
        mEdge makeIdent(QuantumRegisterCount n) {
          return makeIdent(0, static_cast<QuantumRegister>(n - 1));
        }

        mEdge makeIdent(QuantumRegisterCount leastSignificantQubit,
                        QuantumRegisterCount mostSignificantQubit) {
          if (mostSignificantQubit < leastSignificantQubit) {
            return mEdge::one;
          }

          if (leastSignificantQubit == 0 &&
              idTable.at(mostSignificantQubit).nextNode != nullptr) {
            return idTable.at(mostSignificantQubit);
          }

          if (mostSignificantQubit >= 1 &&
              (idTable.at(mostSignificantQubit - 1)).nextNode != nullptr) {
            auto basicDimMost = registersSizes.at(mostSignificantQubit);
            std::vector<mEdge> identityEdges{};

            for (auto i = 0L; i < basicDimMost; i++) {
              for (auto j = 0L; j < basicDimMost; j++) {
                if (i == j) {
                  identityEdges.push_back(idTable[mostSignificantQubit - 1]);
                } else {
                  identityEdges.push_back(mEdge::zero);
                }
              }
            }
            idTable.at(mostSignificantQubit) =
                makeDDNode(static_cast<QuantumRegister>(mostSignificantQubit),
                           identityEdges);

            return idTable.at(mostSignificantQubit);
          }

          // create an Identity DD from scratch
          auto basicDimLeast = registersSizes.at(leastSignificantQubit);
          std::vector<mEdge> identityEdgesLeast{};

          for (auto i = 0L; i < basicDimLeast; i++) {
            for (auto j = 0L; j < basicDimLeast; j++) {
              if (i == j) {
                identityEdgesLeast.push_back(mEdge::one);
              } else {
                identityEdgesLeast.push_back(mEdge::zero);
              }
            }
          }

          auto e =
              makeDDNode(static_cast<QuantumRegister>(leastSignificantQubit),
                         identityEdgesLeast);

          for (std::size_t intermediaryRegs = leastSignificantQubit + 1;
               intermediaryRegs <=
               std::make_unsigned_t<QuantumRegister>(mostSignificantQubit);
               intermediaryRegs++) {
            auto basicDimInt = registersSizes.at(intermediaryRegs);
            std::vector<mEdge> identityEdgesInt{};

            for (auto i = 0L; i < basicDimInt; i++) {
              for (auto j = 0L; j < basicDimInt; j++) {
                if (i == j) {
                  identityEdgesInt.push_back(e);
                } else {
                  identityEdgesInt.push_back(mEdge::zero);
                }
              }
            }
            e = makeDDNode(static_cast<QuantumRegister>(intermediaryRegs),
                           identityEdgesInt);
          }

          if (leastSignificantQubit == 0) {
            idTable.at(mostSignificantQubit) = e;
          }
          return e;
        }

        // identity table access and reset
        [[nodiscard]] const auto& getIdentityTable() const { return idTable; }

        void clearIdentityTable() {
          for (auto& entry : idTable) {
            entry.nextNode = nullptr;
          }
        }

       private:
        std::vector<mEdge> idTable{};

        ///
        /// Vector and matrix extraction from DDs
        ///
       public:
        /// Get a single element of the vector or matrix represented by the dd
        /// with root edge e \tparam Edge type of edge to use (vector or matrix)
        /// \param e edge to traverse
        /// \param path_elements string {0, 1, 2, 3}^n describing which outgoing
        /// edge should be followed
        ///        (for vectors entries are limited to 0 and 1)
        ///        If string is longer than required, the additional characters are ignored.
        /// \return the complex amplitude of the specified element

        template<class Edge>
        ComplexValue getValueByPath(const Edge& edge, const std::string& path_elements) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }

            auto temp_comp_numb = complexNumber.getTemporary(1, 0);
            auto current_edge = edge;
            do {
              ComplexNumbers::mul(temp_comp_numb, temp_comp_numb,
                                  current_edge.weight);
              std::size_t tmp =
                  path_elements.at(current_edge.nextNode->var_indx) - '0';
              assert(tmp <= current_edge.nextNode->edges.size());
              current_edge = current_edge.nextNode->edges.at(tmp);
            } while (!current_edge.isTerminal());

            ComplexNumbers::mul(temp_comp_numb, temp_comp_numb, current_edge.weight);

            return {CTEntry::val(temp_comp_numb.real), CTEntry::val(temp_comp_numb.img)};
        }



        ComplexValue getValueByPath(const vEdge& edge, std::vector<unsigned long>& repr_i) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }
            return getValueByPath(edge, Complex::one, repr_i);
        }



        ComplexValue getValueByPath(const vEdge& edge, const Complex& amp, std::vector<unsigned long>& repr) {
          auto c_numb = complexNumber.mulCached(edge.weight, amp);

          if (edge.isTerminal()) {
            complexNumber.returnToCache(c_numb);
            return {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
          }

          ComplexValue return_amp{};

          if (!edge.nextNode->edges.at(repr.front())
                   .weight.approximatelyZero()) {
            std::vector<unsigned long> reprSlice(repr.begin() + 1, repr.end());
            return_amp = getValueByPath(edge.nextNode->edges.at(repr.front()),
                                        c_numb, reprSlice);
          }

          complexNumber.returnToCache(c_numb);
          return return_amp;
        }

        ComplexValue getValueByPath(const mEdge& edge,
                                    std::vector<unsigned long>& reprI,
                                    std::vector<unsigned long>& reprJ) {
          if (edge.isTerminal()) {
            return {CTEntry::val(edge.weight.real),
                    CTEntry::val(edge.weight.img)};
          }
          return getValueByPath(edge, Complex::one, reprI, reprJ);
        }

        ComplexValue getValueByPath(const mEdge& edge, const Complex& amp,
                                    std::vector<unsigned long>& repr_i,
                                    std::vector<unsigned long>& repr_j) {
          // row major encoding

          auto c_numb = complexNumber.mulCached(edge.weight, amp);

          if (edge.isTerminal()) {
            complexNumber.returnToCache(c_numb);
            return {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
          }

          // for every qubit checks if you're in the correct row or correct
          // column
          // const bool row = i & (1ULL << edge.nextNode->var_indx); // adds one
          // to var indx const bool col = j & (1ULL << edge.nextNode->var_indx);

          const auto row = repr_i.front();
          const auto col = repr_j.front();
          const auto rowMajorIndex = row * edge.nextNode->edges.size() + col;
          ComplexValue returnAmp{};

          if (!edge.nextNode->edges.at(rowMajorIndex)
                   .weight.approximatelyZero()) {
            std::vector<unsigned long> reprSliceI(repr_i.begin() + 1,
                                                  repr_i.end());
            std::vector<unsigned long> reprSliceJ(repr_j.begin() + 1,
                                                  repr_j.end());
            returnAmp = getValueByPath(edge.nextNode->edges.at(rowMajorIndex),
                                       c_numb, reprSliceI, reprSliceJ);
          }
          complexNumber.returnToCache(c_numb);
          return returnAmp;
        }

        CVec getVector(const vEdge& edge) {
          const std::size_t dim = static_cast<const size_t>(
              std::accumulate(registersSizes.begin(), registersSizes.end(), 1,
                              std::multiplies<>()));
          // allocate resulting vector
          auto vec = CVec(dim, {0.0, 0.0});
          getVector(edge, Complex::one, 0, vec, dim);
          return vec;
        }

        void getVector(const vEdge& edge, const Complex& amp, std::size_t i, CVec& vec, std::size_t next) {
            // calculate new accumulated amplitude
            auto cNumb = complexNumber.mulCached(edge.weight, amp);

            // base case
            if (edge.isTerminal()) {
              vec.at(i) = {CTEntry::val(cNumb.real), CTEntry::val(cNumb.img)};
              complexNumber.returnToCache(cNumb);
              return;
            }

            auto offset = (next - i) / edge.nextNode->edges.size();

            for (auto k = 0L; k < edge.nextNode->edges.size(); k++) {
              if (!edge.nextNode->edges.at(k).weight.approximatelyZero()) {
                getVector(edge.nextNode->edges.at(k), cNumb, i + (k * offset),
                          vec, i + ((k + 1) * offset));
              }
            }

            complexNumber.returnToCache(cNumb);
        }

        std::vector<unsigned long> getReprOfIndex(
            const unsigned long i, const unsigned long num_entries) {
          std::vector<unsigned long> repr;
          repr.reserve(numberOfQuantumRegisters);
          // get representation
          auto iIndex = i;
          auto pathWay = 0UL;
          auto cardinality = num_entries;

          auto counter = 0UL;
          auto index = 0UL;

          while (counter < numberOfQuantumRegisters) {
            index = numberOfQuantumRegisters - counter - 1;
            cardinality = cardinality / registersSizes.at(index);
            pathWay = iIndex / cardinality;
            iIndex = iIndex % cardinality;

            repr.push_back(pathWay);
            counter = counter + 1;
          }

          return repr;
        }

        void printVector(const vEdge& edge) {
          unsigned long long numEntries = static_cast<unsigned long long int>(
              std::accumulate(registersSizes.begin(), registersSizes.end(), 1,
                              std::multiplies<>()));

          for (auto i = 0ULL; i < numEntries; i++) {
            auto reprI = getReprOfIndex(i, numEntries);
            // get amplitude
            const auto amplitude = getValueByPath(edge, reprI);
            // TODO HOW SHALL WE REPRESENT??
            //
            for (const unsigned long& coeff : reprI) {
              std::cout << coeff;
            }
            reprI.clear();

            constexpr auto precision = 3;
            // set fixed width to maximum of a printed number
            // (-) 0.precision plus/minus 0.precision i
            constexpr auto width = 1 + 2 + precision + 1 + 2 + precision + 1;
            std::cout << ": " << std::setw(width)
                      << ComplexValue::toString(amplitude.r, amplitude.i, false,
                                                precision)
                      << "\n";
          }
          std::cout << std::flush;
        }


    private:

        // check whether node represents a symmetric matrix or the identity
        void checkSpecialMatrices(mNode* node) {
       if (node->varIndx == -1) return;

       node->identity = false;   // assume not identity
       node->symmetric = false;  // assume symmetric

       // check if matrix is symmetric
       auto numberOfEdges = node->edges.size();

       auto basicDim = registersSizes.at(node->varIndx);

       for (auto i = 0L; i < basicDim; i++) {
         if (!node->edges.at(i * basicDim + i).nextNode->symmetric) {
           return;
         }
       }
       // TODO WHY RETURN IF DIAGONAL IS SYMMETRIC??
       // if (!node->edges.at(0).nextNode->symmetric ||
       // !node->edges.at(3).nextNode->symmetric) return;

       for (auto i = 0L; i < basicDim; i++) {
         for (auto j = 0L; j < basicDim; j++) {
           if (i != j) {
             // row major indexing - enable optimization here
             if (transpose(node->edges.at(i * basicDim + j)) !=
                 node->edges.at(j * basicDim + i)) {
               return;
             }
           }
         }
       }
       // if (transpose(node->edges.at(1)) != node->edges.at(2)) return;

       node->symmetric = true;

       // check if matrix resembles identity
       for (auto i = 0L; i < basicDim; i++) {
         for (auto j = 0L; j < basicDim; j++) {
           // row major indexing - enable optimization here
           if (i == j) {
             if (!(node->edges[i * basicDim + j].nextNode->identity) ||
                 (node->edges[i * basicDim + j].weight) != Complex::one)
               return;
           } else {
             if ((node->edges[i * basicDim + j].weight) != Complex::zero)
               return;
           }
         }
       }
       /*
       if (!(p->e[0].p->ident) || (p->e[1].w) != Complex::zero ||
           (p->e[2].w) != Complex::zero || (p->e[0].w) != Complex::one ||
           (p->e[3].w) != Complex::one || !(p->e[3].p->ident))
           return;
       */
       node->identity = true;
     }

     ///
        /// Matrix (conjugate) transpose
        ///
    public:
        //todo figure out the parameters here
        UnaryComputeTable<mEdge, mEdge, 4096> matrixTranspose{};
        UnaryComputeTable<mEdge, mEdge, 4096> conjugateMatrixTranspose{};

        mEdge transpose(const mEdge& edge) {
          if (edge.nextNode == nullptr || edge.isTerminal() ||
              edge.nextNode->symmetric) {
            return edge;
          }

          // check in compute table
          auto result = matrixTranspose.lookup(edge);
          if (result.nextNode != nullptr) {
            return result;
          }

          std::vector<mEdge> newEdge{};
          auto basicDim = registersSizes.at(edge.nextNode->varIndx);

          // transpose sub-matrices and rearrange as required
          for (auto i = 0U; i < basicDim; i++) {
            for (auto j = 0U; j < basicDim; j++) {
              newEdge.at(basicDim * i + j) =
                  transpose(edge.nextNode->edges.at(basicDim * j + i));
            }
          }
          // create new top node
          result = makeDDNode(edge.nextNode->varIndx, newEdge);
          // adjust top weight
          auto c = complexNumber.getTemporary();
          ComplexNumbers::mul(c, result.weight, edge.weight);
          result.weight = complexNumber.lookup(c);

          // put in compute table
          matrixTranspose.insert(edge, result);
          return result;
        }
        mEdge conjugateTranspose(const mEdge& edge) {
          if (edge.nextNode == nullptr) return edge;
          if (edge.isTerminal()) {  // terminal case
            auto result = edge;
            result.weight = ComplexNumbers::conj(edge.weight);
            return result;
          }

          // check if in compute table
          auto result = conjugateMatrixTranspose.lookup(edge);
          if (result.nextNode != nullptr) {
            return result;
          }

          std::vector<mEdge> newEdge{};
          auto basicDim = registersSizes.at(edge.nextNode->varIndx);

          // conjugate transpose submatrices and rearrange as required
          for (auto i = 0U; i < basicDim; ++i) {
            for (auto j = 0U; j < basicDim; ++j) {
              newEdge.at(basicDim * i + j) =
                  conjugateTranspose(edge.nextNode->edges.at(basicDim * j + i));
            }
          }
          // create new top node
          result = makeDDNode(edge.nextNode->varIndx, newEdge);

          auto c = complexNumber.getTemporary();
          // adjust top weight including conjugate
          ComplexNumbers::mul(c, result.weight,
                              ComplexNumbers::conj(edge.weight));
          result.weight = complexNumber.lookup(c);

          // put it in the compute table
          conjugateMatrixTranspose.insert(edge, result);
          return result;
        }

        ///
        /// Unique tables, Reference counting and garbage collection
        ///
    public:
        // unique tables
        template<class Node>
        [[nodiscard]] UniqueTable<Node> &getUniqueTable();

     template <class Node>
     void incRef(const Edge<Node>& e) {
       getUniqueTable<Node>().incRef(e);
     }

     template <class Node>
     void decRef(const Edge<Node>& e) {
       getUniqueTable<Node>().decRef(e);
     }

     UniqueTable<vNode> vUniqueTable{numberOfQuantumRegisters};
     UniqueTable<mNode> mUniqueTable{numberOfQuantumRegisters};
    };

    void clearUniqueTables() {
        // TODO IMPLEMENT
        //vUniqueTable.clear();
        //mUniqueTable.clear();
    }


    inline MDDPackage::vNode MDDPackage::vNode::terminalNode{
            {
                    {
                            {
                                    nullptr, Complex::zero
                            }, {
                                    nullptr, Complex::zero
                            }
                    }
            },
            nullptr,
            0,
            -1
    };

    inline MDDPackage::mNode MDDPackage::mNode::terminalNode{
            {
                    {
                            {
                                    nullptr, Complex::zero
                            }, {
                                    nullptr, Complex::zero
                            }, {
                                    nullptr, Complex::zero
                            }, {
                                    nullptr, Complex::zero
                            }
                    }
            },
            nullptr,
            0,
            -1,
            true,
            true};

    template<>
    [[nodiscard]] inline UniqueTable<MDDPackage::vNode> &MDDPackage::getUniqueTable() { return vUniqueTable; }

    template<>
    [[nodiscard]] inline UniqueTable<MDDPackage::mNode> &MDDPackage::getUniqueTable() { return mUniqueTable; }
    /*
    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::vCachedEdge,
    MDDPackage::vCachedEdge, MDDPackage::vCachedEdge>&
    MDDPackage::getAddComputeTable() { return vectorAdd; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::mCachedEdge,
    MDDPackage::mCachedEdge, MDDPackage::mCachedEdge>&
    MDDPackage::getAddComputeTable() { return matrixAdd; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::mEdge, MDDPackage::vEdge,
    MDDPackage::vCachedEdge>& MDDPackage::getMultiplicationComputeTable() {
    return matrixVectorMultiplication; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::mEdge, MDDPackage::mEdge,
    MDDPackage::mCachedEdge>& MDDPackage::getMultiplicationComputeTable() {
    return matrixMatrixMultiplication; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::vEdge, MDDPackage::vEdge,
    MDDPackage::vCachedEdge, 4096>& MDDPackage::getKroneckerComputeTable() {
    return vectorKronecker; }

    template<>
    [[nodiscard]] inline ComputeTable<Package::mEdge, Package::mEdge,
    Package::mCachedEdge, 4096>& Package::getKroneckerComputeTable() { return
    matrixKronecker; }
    */
    }  // namespace dd

#endif