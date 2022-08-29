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


#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
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
        ComplexNumbers complex_number{};

        ///
        /// Construction, destruction, information and reset
        ///
    public:


        explicit MDDPackage(std::size_t nqr, std::vector<size_t> sizes):
        complex_number(ComplexNumbers()), number_of_quantum_registers(nqr), registers_sizes(std::move(sizes)) {};

        ~MDDPackage()                      = default;
        MDDPackage(const MDDPackage& MDDPackage) = delete; // no copy constructor
        MDDPackage& operator=(const MDDPackage& MDDPackage) = delete; // no copy assignment constructor

        //TODO RESIZE

        //TODO CHECK SYNTAX OF SETTERS AND GETTERS
        // getter for number qudits

        [[nodiscard]] auto qregisters() const { return number_of_quantum_registers; }

        // setter dimensionalisties
        [[nodiscard]] auto register_dimensions( const std::vector<size_t>& regs ) {
            registers_sizes =  regs;
        }

        // getter for sizes
        [[nodiscard]] auto regs_size() const { return registers_sizes; }

    private:
        std::size_t number_of_quantum_registers;
        //TODO THIS IS NOT CONST RIGHT?
        std::vector<size_t> registers_sizes;

        ///
        /// Vector nodes, edges and quantum states
        ///
    public:
        struct vNode {
            std::vector< Edge<vNode> >     edges_{};    // edges out of this node
            vNode*                         next_{}; // used to link nodes in unique table
            RefCount                       ref_count{};  // reference count
            QuantumRegister                var_indx{};    // variable index (nonterminal) value (-1 for terminal)

            static vNode            terminalNode;
            constexpr static vNode* terminal{&terminalNode};

            static constexpr bool isTerminal(const vNode* node_point) { return node_point == terminal; }
        };
        using vEdge       = Edge<vNode>;
        using vCachedEdge = CachedEdge<vNode>;





        vEdge normalize(const vEdge& edge, bool cached) {

            std::vector< bool > zero;

            for (auto & i : edge.next_node->edges_){
                zero.push_back(i.weight.approximatelyZero());
            }

            // make sure to release cached numbers approximately zero, but not exactly zero
            if (cached) {
                for (auto i = 0U; i < zero.size(); i++) {

                    if (zero.at(i) && edge.next_node->edges_.at(i).weight != Complex::zero) {
                        //TODO what is returnToCache

                        complex_number.returnToCache(edge.next_node->edges_.at(i).weight);
                        edge.next_node->edges_.at(i) = vEdge::zero;
                    }
                }
            }

            // all equal to zero
            if ( none_of( cbegin(zero), cend(zero), std::logical_not<bool>() ) ) {

                if (!cached && !edge.isTerminal()) {
                    // If it is not a cached computation, the node has to be put back into the chain
                    vUniqueTable.returnNode(edge.next_node);
                }
                return vEdge::zero;
            }

            // find indices that are not zero
            std::vector<int> non_zero_indices;
            for(auto i = 0U; i < zero.size(); i++) {
                if(!zero.at(i)) {
                    non_zero_indices.push_back(i);
                }
            }

            if ( non_zero_indices.size() == 1) {
                // search for first element different from zero
                auto  current_edge = edge;
                auto& weight_from_child = current_edge.next_node->edges_[non_zero_indices.front()].weight;

                if(cached && weight_from_child != Complex::one) {
                    current_edge.weight = weight_from_child;
                } else {
                    current_edge.weight = complex_number.lookup(weight_from_child);
                }

                weight_from_child = Complex::one;
                return current_edge;
            }

            //calculate normalizing factor
            auto sum_norm2 = ComplexNumbers::mag2(edge.next_node->edges_.at(0).weight);

            for(auto i = 1; i < edge.next_node->edges_.size(); i++){
                sum_norm2 = sum_norm2 + ComplexNumbers::mag2(edge.next_node->edges_.at(i).weight);
            }

            const auto commonFactor         = std::sqrt(sum_norm2);



            // set incoming edge weight to max
            auto  current_edge   = edge;


            //TODO CHECK EXACTLY THIS CACHING SYSTEM
            //auto& max = r.next_node->edges_[argMax];
            if (cached && current_edge.weight != Complex::one) {
                //current_edge.weight = max.weight;
                current_edge.weight.real->value *= commonFactor;
                current_edge.weight.img->value *= commonFactor;
            } else {
                current_edge.weight = complex_number.lookup(CTEntry::val(current_edge.weight.real) * commonFactor, CTEntry::val(current_edge.weight.img) * commonFactor);
                if (current_edge.weight.approximatelyZero()) {
                    return vEdge::zero;
                }
            }

            // actual normalization of the edges
            for (auto & i : edge.next_node->edges_){
                if (cached) {
                    complex_number.returnToCache(i.weight);
                    ComplexNumbers::div(i.weight, i.weight, current_edge.weight);
                    i.weight = complex_number.lookup(i.weight);
                } else {
                    auto c = complex_number.getTemporary();
                    ComplexNumbers::div(c, i.weight, current_edge.weight);
                    i.weight = complex_number.lookup(c);
                }
                if (i.weight == Complex::zero) {
                    i = vEdge::zero;
                }
            }

            return current_edge;
        }


        // generate |0...0> with N quantum registers
        vEdge makeZeroState(QuantumRegisterCount n, std::size_t start = 0) {
            if (n + start > number_of_quantum_registers) {
                //TODO UNDERSTAND RESIZING
                throw std::runtime_error("Requested state with " +
                                         std::to_string(n + start) +
                                         " QUANTUM REGISTERS, but current package configuration only supports up to " +
                                         std::to_string(number_of_quantum_registers) +
                                         " QUANTUM REGISTERS. Please allocate a larger package instance.");
            }
            auto first = vEdge::one;
            for (std::size_t node_idx = start; node_idx < n + start; node_idx++) {

                std::vector<Edge<vNode>> new_outgoing_edges ;
                new_outgoing_edges.reserve(registers_sizes.at(node_idx));
                new_outgoing_edges.push_back(first);
                for(auto i = 1U; i < registers_sizes.at(node_idx); i++){
                    new_outgoing_edges.push_back(vEdge::zero);
                }

                first = makeDDNode(static_cast<QuantumRegister>(node_idx), new_outgoing_edges );
            }
            return first;
        }


        // create a normalized DD node and return an edge pointing to it. The node is not recreated if it already exists.
        template<class Node>
        Edge<Node> makeDDNode(QuantumRegister varidx, const std::vector<Edge<Node>>& edges, bool cached = false) {

            auto&      uniqueTable = getUniqueTable<Node>();

            Edge<Node> new_edge { uniqueTable.getNode(), Complex::one };
            new_edge.next_node->var_indx = varidx;
            new_edge.next_node->edges_ = edges;

            assert(new_edge.next_node->ref_count == 0);

            for ([[maybe_unused]] const auto& edge: edges)
                assert( edge.next_node->var_indx == varidx - 1 || edge.isTerminal() );

            // normalize it
            new_edge = normalize(new_edge, cached);
            assert(new_edge.next_node->var_indx == varidx || new_edge.isTerminal());

            // look it up in the unique tables
            auto looked_up_edge = uniqueTable.lookup(new_edge, false);
            assert(looked_up_edge.next_node->var_indx == varidx || looked_up_edge.isTerminal());

            // set specific node properties for matrices
            if constexpr (std::is_same_v<Node, mNode>){
                if (looked_up_edge.next_node == new_edge.next_node)
                    checkSpecialMatrices(looked_up_edge.next_node);
            }

            return looked_up_edge;
        }







    public:
        struct mNode {
            std::vector< Edge<mNode> >     edges_{};           // edges out of this node
            mNode*                         next_{};        // used to link nodes in unique table
            RefCount                       ref_count{};         // reference count
            QuantumRegister                var_indx{};           // variable index (nonterminal) value (-1 for terminal)
            bool                           symmetric  = false; // node is symmetric
            bool                           identity = false; // node resembles identity

            static mNode            terminalNode;
            constexpr static mNode* terminal{&terminalNode};

            static constexpr bool isTerminal(const mNode* node_point) { return node_point == terminal; }
        };
        using mEdge       = Edge<mNode>;
        using mCachedEdge = CachedEdge<mNode>;


        mEdge normalize(const mEdge& edge, bool cached) {
            auto argmax = -1;

            std::vector< bool > zero;

            for (auto & i : edge.next_node->edges_){
                zero.push_back(i.weight.approximatelyZero());
            }

            // make sure to release cached numbers approximately zero, but not exactly zero
            if (cached) {
                for (auto i = 0U; i < zero.size(); i++) {

                    if (zero.at(i) && edge.next_node->edges_.at(i).weight != Complex::zero) {
                        //TODO what is returnToCache

                        complex_number.returnToCache(edge.next_node->edges_.at(i).weight);
                        edge.next_node->edges_.at(i) = mEdge::zero;
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
                    max_magnitude    = ComplexNumbers::mag2(edge.next_node->edges_.at(i).weight);
                    max_weight   = edge.next_node->edges_.at(i).weight;
                } else {
                    auto current_magnitude = ComplexNumbers::mag2(edge.next_node->edges_.at(i).weight);
                    if (current_magnitude - max_magnitude > ComplexTable<>::tolerance()) {
                        argmax = static_cast<decltype(argmax)>(i);
                        max_magnitude    = current_magnitude;
                        max_weight   = edge.next_node->edges_.at(i).weight;
                    }
                }
            }

            // all equal to zero
            if (argmax == -1) {
                if (!cached && !edge.isTerminal()) {
                    // If it is not a cached computation, the node has to be put back into the chain
                    mUniqueTable.returnNode(edge.next_node);
                }
                return mEdge::zero;
            }

            auto current_edge = edge;
            // divide each entry by max
            for (auto i = 0U; i < edge.next_node->edges_.size(); ++i) {
                if (static_cast<decltype(argmax)>(i) == argmax) {
                    if (cached) {
                        if (current_edge.weight == Complex::one)
                            current_edge.weight = max_weight;
                        else
                            ComplexNumbers::mul(current_edge.weight, current_edge.weight, max_weight);
                    } else {
                        if (current_edge.weight == Complex::one) {
                            current_edge.weight = max_weight;
                        } else {
                            auto new_complex_numb = complex_number.getTemporary();
                            ComplexNumbers::mul(new_complex_numb, current_edge.weight, max_weight);
                            current_edge.weight = complex_number.lookup(new_complex_numb);
                        }
                    }
                    current_edge.next_node->edges_.at(i).weight = Complex::one;
                } else {
                    if (cached && !zero.at(i) && current_edge.next_node->edges_.at(i).weight != Complex::one) {
                        complex_number.returnToCache(current_edge.next_node->edges_.at(i).weight);
                    }
                    if (current_edge.next_node->edges_.at(i).weight.approximatelyOne())
                        current_edge.next_node->edges_.at(i).weight = Complex::one;
                    auto new_complex_numb = complex_number.getTemporary();

                    ComplexNumbers::div(new_complex_numb, current_edge.next_node->edges_.at(i).weight, max_weight);
                    current_edge.next_node->edges_.at(i).weight = complex_number.lookup(new_complex_numb);
                }
            }
            return current_edge;
        }

        ///
        /// Vector and matrix extraction from DDs
        ///
    public:
        /// Get a single element of the vector or matrix represented by the dd with root edge e
        /// \tparam Edge type of edge to use (vector or matrix)
        /// \param e edge to traverse
        /// \param path_elements string {0, 1, 2, 3}^n describing which outgoing edge should be followed
        ///        (for vectors entries are limited to 0 and 1)
        ///        If string is longer than required, the additional characters are ignored.
        /// \return the complex amplitude of the specified element

        template<class Edge>
        ComplexValue getValueByPath(const Edge& edge, const std::string& path_elements) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }

            auto temp_comp_numb = complex_number.getTemporary(1, 0);
            auto current_edge = edge;
            do {
                ComplexNumbers::mul(temp_comp_numb, temp_comp_numb, current_edge.weight);
                std::size_t tmp = path_elements.at(current_edge.next_node->var_indx) - '0';
                assert(tmp <= current_edge.next_node->edges_.size());
                current_edge = current_edge.next_node->edges_.at(tmp);
            } while (!current_edge.isTerminal());

            ComplexNumbers::mul(temp_comp_numb, temp_comp_numb, current_edge.weight);

            return {CTEntry::val(temp_comp_numb.real), CTEntry::val(temp_comp_numb.img)};
        }



        ComplexValue getValueByPath(const vEdge& edge, std::size_t i) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }
            return getValueByPath(edge, Complex::one, i);
        }



        ComplexValue getValueByPath(const vEdge& edge, const Complex& amp, std::size_t i) {
            auto c_numb = complex_number.mulCached(edge.weight, amp);

            if (edge.isTerminal()) {
                complex_number.returnToCache(c_numb);
                return {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
            }

            const bool one = i & (1ULL << edge.next_node->var_indx);

            ComplexValue return_amp{};

            if (!one && !edge.next_node->edges_[0].weight.approximatelyZero()) {
                return_amp = getValueByPath(edge.next_node->edges_[0], c_numb, i);
            } else if (one && !edge.next_node->edges_[1].weight.approximatelyZero()) {
                return_amp = getValueByPath(edge.next_node->edges_[1], c_numb, i);
            }
            complex_number.returnToCache(c_numb);
            return return_amp;
        }


        ComplexValue getValueByPath(const mEdge& edge, std::size_t i, std::size_t j) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }
            return getValueByPath(edge, Complex::one, i, j);
        }



        ComplexValue getValueByPath(const mEdge& edge, const Complex& amp, std::size_t i, std::size_t j) {
            auto c_numb = complex_number.mulCached(edge.weight, amp);

            if (edge.isTerminal()) {
                complex_number.returnToCache(c_numb);
                return {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
            }

            const bool row = i & (1ULL << edge.next_node->var_indx);
            const bool col = j & (1ULL << edge.next_node->var_indx);

            ComplexValue r{};

            if (!row && !col && !edge.next_node->edges_[0].weight.approximatelyZero()) {
                r = getValueByPath(edge.next_node->edges_[0], c_numb, i, j);
            } else if (!row && col && !edge.next_node->edges_[1].weight.approximatelyZero()) {
                r = getValueByPath(edge.next_node->edges_[1], c_numb, i, j);
            } else if (row && !col && !edge.next_node->edges_[2].weight.approximatelyZero()) {
                r = getValueByPath(edge.next_node->edges_[2], c_numb, i, j);
            } else if (row && col && !edge.next_node->edges_[3].weight.approximatelyZero()) {
                r = getValueByPath(edge.next_node->edges_[3], c_numb, i, j);
            }
            complex_number.returnToCache(c_numb);
            return r;
        }



        CVec getVector(const vEdge& edge) {
            const std::size_t dim = 2ULL << edge.next_node->var_indx;
            // allocate resulting vector
            auto vec = CVec(dim, {0.0, 0.0});
            getVector(edge, Complex::one, 0, vec);
            return vec;
        }


        void getVector(const vEdge& edge, const Complex& amp, std::size_t i, CVec& vec) {
            // calculate new accumulated amplitude
            auto c_numb = complex_number.mulCached(edge.weight, amp);

            // base case
            if (edge.isTerminal()) {
                vec.at(i) = {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
                complex_number.returnToCache(c_numb);
                return;
            }

            const std::size_t x = i | (1ULL << edge.next_node->var_indx);

            // recursive case
            if (!edge.next_node->edges_[0].weight.approximatelyZero())
                getVector(edge.next_node->edges_[0], c_numb, i, vec);
            if (!edge.next_node->edges_[1].weight.approximatelyZero())
                getVector(edge.next_node->edges_[1], c_numb, x, vec);
            complex_number.returnToCache(c_numb);
        }


        //TODO WHY THERE ARE COMPLEX VALUES COMPLES NUMBERS AND JUST COMPLEX

        void printVector(const vEdge& edge) {
            unsigned long long length_vec = std::accumulate(registers_sizes.begin(),registers_sizes.end(), 1, std::multiplies<>());
            const unsigned long long element = length_vec << edge.next_node->var_indx;

            for (auto i = 0ULL; i < element; i++) {
                const auto amplitude = getValueByPath(edge, i);

                for (QuantumRegister j = edge.next_node->var_indx; j >= 0; j--) {
                    std::cout << ((i >> j) & 1ULL);
                }
                constexpr auto precision = 3;
                // set fixed width to maximum of a printed number
                // (-) 0.precision plus/minus 0.precision i
                constexpr auto width = 1 + 2 + precision + 1 + 2 + precision + 1;
                std::cout << ": " << std::setw(width) << ComplexValue::toString(amplitude.r, amplitude.i, false, precision) << "\n";
            }
            std::cout << std::flush;
        }



        ///
        /// Unique tables, Reference counting and garbage collection
        ///
    public:
        // unique tables
        template<class Node>
        [[nodiscard]] UniqueTable<Node>& getUniqueTable();

        template<class Node>
        void incRef(const Edge<Node>& e) {
            getUniqueTable<Node>().incRef(e);
        }
        template<class Node>
        void decRef(const Edge<Node>& e) {
            getUniqueTable<Node>().decRef(e);
        }

        UniqueTable<vNode> vUniqueTable{number_of_quantum_registers};
        UniqueTable<mNode> mUniqueTable{number_of_quantum_registers};
    };

    inline MDDPackage::vNode MDDPackage::vNode::terminalNode{{{{nullptr, Complex::zero}, {nullptr, Complex::zero}}},
                                                       nullptr,
                                                       0,
                                                       -1};

    inline MDDPackage::mNode MDDPackage::mNode::terminalNode{
            {{{nullptr, Complex::zero}, {nullptr, Complex::zero}, {nullptr, Complex::zero}, {nullptr, Complex::zero}}},
            nullptr,
            0,
            -1,
            true,
            true};

    template<>
    [[nodiscard]] inline UniqueTable<MDDPackage::vNode>& MDDPackage::getUniqueTable() { return vUniqueTable; }

    template<>
    [[nodiscard]] inline UniqueTable<MDDPackage::mNode>& MDDPackage::getUniqueTable() { return mUniqueTable; }
    /*
    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::vCachedEdge, MDDPackage::vCachedEdge, MDDPackage::vCachedEdge>& MDDPackage::getAddComputeTable() { return vectorAdd; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::mCachedEdge, MDDPackage::mCachedEdge, MDDPackage::mCachedEdge>& MDDPackage::getAddComputeTable() { return matrixAdd; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::mEdge, MDDPackage::vEdge, MDDPackage::vCachedEdge>& MDDPackage::getMultiplicationComputeTable() { return matrixVectorMultiplication; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::mEdge, MDDPackage::mEdge, MDDPackage::mCachedEdge>& MDDPackage::getMultiplicationComputeTable() { return matrixMatrixMultiplication; }

    template<>
    [[nodiscard]] inline ComputeTable<MDDPackage::vEdge, MDDPackage::vEdge, MDDPackage::vCachedEdge, 4096>& MDDPackage::getKroneckerComputeTable() { return vectorKronecker; }

    template<>
    [[nodiscard]] inline ComputeTable<Package::mEdge, Package::mEdge, Package::mCachedEdge, 4096>& Package::getKroneckerComputeTable() { return matrixKronecker; }
    */
}


#endif