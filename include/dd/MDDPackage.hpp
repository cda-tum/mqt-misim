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

        ~MDDPackage()                      = default;
        MDDPackage(const MDDPackage& MDDPackage) = delete; // no copy constructor
        MDDPackage& operator=(const MDDPackage& MDDPackage) = delete; // no copy assignment constructor


    private:
        std::size_t number_of_quantum_registers;

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
            if ( none_of(cbegin(zero), cend(zero), std::logical_not<bool>()) ) {

                if (!cached && !edge.isTerminal()) {
                    // If it is not a cached computation, the node has to be put back into the chain
                    vUniqueTable.returnNode(edge.next_node);
                }
                return vEdge::zero;
            }

            // exactly one is non zero
            auto nonZero = 0U;
            auto nZindex = 0U;
            for(auto i = 0U; i < zero.size(); i++) {
                if (!zero.at(i)) {
                    nonZero++;
                    nZindex = i;
                }
            }
            if ( nonZero == 1) {
                // search for first element different from zero
                auto  r = edge;
                auto& w = r.next_node->edges_[nZindex].weight;
                if (cached && w != Complex::one) {
                    r.weight = w;
                } else {
                    r.weight = complex_number.lookup(w);
                }
                w = Complex::one;
                return r;
            }

            //calculate normalizing factor
            auto norm2 = ComplexNumbers::mag2(edge.next_node->edges_.at(0).weight);

            for(auto i = 1U; i < edge.next_node->edges_.size(); i++){
                norm2 = norm2 + ComplexNumbers::mag2(edge.next_node->edges_.at(i).weight);
            }
            //const auto norm2        = mag0 + mag1;
            auto mag2Max = ComplexNumbers::mag2(edge.next_node->edges_.at(1).weight);
            auto argMax  = 1U;
            for(auto i = 0U; i < edge.next_node->edges_.size(); i++){
                auto mag0 = ComplexNumbers::mag2(edge.next_node->edges_.at(i-1).weight);

                mag2Max     = ( mag0 + ComplexTable<>::tolerance() >= mag2Max  ) ? mag0  : mag2Max;
                argMax      = ( mag0 + ComplexTable<>::tolerance() >= mag2Max  ) ? i  : argMax;
            }
            //const auto mag2Max      = (mag0 + ComplexTable<>::tolerance() >= mag1) ? mag0 : mag1;
            //const auto argMax       = (mag0 + ComplexTable<>::tolerance() >= mag1) ? 0 : 1;
            const auto norm         = std::sqrt(norm2);
            const auto magMax       = std::sqrt(mag2Max);
            const auto commonFactor = norm / magMax;


            // set incoming edge weight to max
            auto  r   = edge;
            auto& max = r.next_node->edges_[argMax];
            if (cached && max.weight != Complex::one) {
                r.weight = max.weight;
                r.weight.real->value *= commonFactor;
                r.weight.img->value *= commonFactor;
            } else {
                r.weight = complex_number.lookup(CTEntry::val(max.weight.real) * commonFactor, CTEntry::val(max.weight.img) * commonFactor);
                if (r.weight.approximatelyZero()) {
                    return vEdge::zero;
                }
            }
            //

            max.weight = complex_number.lookup(magMax / norm, 0.);
            if (max.weight == Complex::zero)
                max = vEdge::zero;

            for (auto & i : edge.next_node->edges_){
                if (cached) {
                    complex_number.returnToCache(i.weight);
                    ComplexNumbers::div(i.weight, i.weight, r.weight);
                    i.weight = complex_number.lookup(i.weight);
                } else {
                    auto c = complex_number.getTemporary();
                    ComplexNumbers::div(c, i.weight, r.weight);
                    i.weight = complex_number.lookup(c);
                }
                if (i.weight == Complex::zero) {
                    i = vEdge::zero;
                }
            }

            return r;
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

        /*
        mEdge normalize(const mEdge& e, bool cached) {
            auto argmax = -1;

            std::vector< bool > zero;

            for (auto & i : e.p->e){
                zero.push_back(i.w.approximatelyZero());
            }

            // make sure to release cached numbers approximately zero, but not exactly zero
            if (cached) {
                for (auto i = 0U; i < e.p->e.size(); i++) {
                    if (zero.at(i) && e.p->e.at(i).w != Complex::zero) {
                        cn.returnToCache(e.p->e.at(i).w);
                        e.p->e.at(i) = mEdge::zero;
                    }
                }
            }

            fp   max  = 0;
            auto maxc = Complex::one;
            // determine max amplitude
            for (auto i = 0U; i < zero.size(); ++i) {
                if (zero.at(i)) continue;
                if (argmax == -1) {
                    argmax = static_cast<decltype(argmax)>(i);
                    max    = ComplexNumbers::mag2(e.p->e.at(i).w);
                    maxc   = e.p->e.at(i).w;
                } else {
                    auto mag = ComplexNumbers::mag2(e.p->e.at(i).w);
                    if (mag - max > ComplexTable<>::tolerance()) {
                        argmax = static_cast<decltype(argmax)>(i);
                        max    = mag;
                        maxc   = e.p->e.at(i).w;
                    }
                }
            }

            // all equal to zero
            if (argmax == -1) {
                if (!cached && !e.isTerminal()) {
                    // If it is not a cached computation, the node has to be put back into the chain
                    mUniqueTable.returnNode(e.p);
                }
                return mEdge::zero;
            }

            auto r = e;
            // divide each entry by max
            for (auto i = 0U; i < e.p->e.size(); ++i) {
                if (static_cast<decltype(argmax)>(i) == argmax) {
                    if (cached) {
                        if (r.w == Complex::one)
                            r.w = maxc;
                        else
                            ComplexNumbers::mul(r.w, r.w, maxc);
                    } else {
                        if (r.w == Complex::one) {
                            r.w = maxc;
                        } else {
                            auto c = cn.getTemporary();
                            ComplexNumbers::mul(c, r.w, maxc);
                            r.w = cn.lookup(c);
                        }
                    }
                    r.p->e.at(i).w = Complex::one;
                } else {
                    if (cached && !zero.at(i) && r.p->e.at(i).w != Complex::one) {
                        cn.returnToCache(r.p->e.at(i).w);
                    }
                    if (r.p->e.at(i).w.approximatelyOne())
                        r.p->e.at(i).w = Complex::one;
                    auto c = cn.getTemporary();

                    ComplexNumbers::div(c, r.p->e.at(i).w, maxc);
                    r.p->e.at(i).w = cn.lookup(c);
                }
            }
            return r;
        }*/
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