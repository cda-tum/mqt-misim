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
            RefCount                       ref_count{};  // reference count, how many active dd are using the node
            QuantumRegister                var_indx{};    // variable index (nonterminal) value (-1 for terminal), index in the circuit endianess 0 from below

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
                if (looked_up_edge.next_node == new_edge.next_node){}

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
        /// Identity matrices
        ///
    public:
        // create n-qudit identity DD. makeIdent(n) === makeIdent(0, n-1)
        mEdge makeIdent(QuantumRegisterCount n) { return makeIdent(0, static_cast<QuantumRegister>(n - 1)); }

        mEdge makeIdent(QuantumRegisterCount leastSignificantQubit, QuantumRegisterCount mostSignificantQubit) {
            if (mostSignificantQubit < leastSignificantQubit)
                return mEdge::one;

            if (leastSignificantQubit == 0 && IdTable.at(mostSignificantQubit).next_node != nullptr) {
                return IdTable.at(mostSignificantQubit);
            }

            if (mostSignificantQubit >= 1 && (IdTable.at(mostSignificantQubit - 1) ).next_node != nullptr) {

                IdTable.at(mostSignificantQubit) = makeDDNode(mostSignificantQubit,
                                                           std::vector{IdTable[mostSignificantQubit - 1],
                                                                      mEdge::zero,
                                                                      mEdge::zero,
                                                                      IdTable[mostSignificantQubit - 1]});
                return IdTable.at(mostSignificantQubit);
            }

            auto e = makeDDNode(leastSignificantQubit, std::vector{mEdge::one, mEdge::zero, mEdge::zero, mEdge::one});

            for (std::size_t k = leastSignificantQubit + 1; k <= std::make_unsigned_t<QuantumRegister>(mostSignificantQubit); k++) {

                e = makeDDNode(static_cast<QuantumRegister>(k), std::vector{e, mEdge::zero, mEdge::zero, e});
            }

            if (leastSignificantQubit == 0)
                IdTable.at(mostSignificantQubit) = e;
            return e;
        }

        // identity table access and reset
        [[nodiscard]] const auto& getIdentityTable() const { return IdTable; }

        void clearIdentityTable() {
            for (auto& entry: IdTable) entry.next_node = nullptr;
        }

    private:
        std::vector<mEdge> IdTable{};

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



        ComplexValue getValueByPath(const vEdge& edge, std::vector<unsigned long>& repr_i) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }
            return getValueByPath(edge, Complex::one, repr_i);
        }



        ComplexValue getValueByPath(const vEdge& edge, const Complex& amp, std::vector<unsigned long>& repr) {
            auto c_numb = complex_number.mulCached(edge.weight, amp);


            if (edge.isTerminal()) {
                complex_number.returnToCache(c_numb);
                return {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
            }

            ComplexValue return_amp{};


            if (!edge.next_node->edges_.at(repr.front()).weight.approximatelyZero()) {
                std::vector<unsigned long> repr_slice(repr.begin()+1, repr.end());
                return_amp = getValueByPath(edge.next_node->edges_.at(repr.front()), c_numb, repr_slice);

            }

            complex_number.returnToCache(c_numb);
            return return_amp;
        }


        ComplexValue getValueByPath(const mEdge& edge, std::vector<unsigned long>& repr_i, std::vector<unsigned long>& repr_j) {
            if (edge.isTerminal()) {
                return {CTEntry::val(edge.weight.real), CTEntry::val(edge.weight.img)};
            }
            return getValueByPath(edge, Complex::one, repr_i, repr_j);
        }



        ComplexValue getValueByPath(const mEdge& edge, const Complex& amp, std::vector<unsigned long>&  repr_i, std::vector<unsigned long>& repr_j) {
            //row major encoding

            auto c_numb = complex_number.mulCached(edge.weight, amp);

            if (edge.isTerminal()) {
                complex_number.returnToCache(c_numb);
                return {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
            }

            // for every qubit checks if you're in the correct row or correct column
            //const bool row = i & (1ULL << edge.next_node->var_indx); // adds one to var indx
            //const bool col = j & (1ULL << edge.next_node->var_indx);

            const auto row = repr_i.front();
            const auto col = repr_j.front();
            const auto row_major_index = row*edge.next_node->edges_.size() + col;
            ComplexValue return_amp{};

            if (!edge.next_node->edges_.at(row_major_index).weight.approximatelyZero()) {
                std::vector<unsigned long> repr_slice_i(repr_i.begin()+1, repr_i.end());
                std::vector<unsigned long> repr_slice_j(repr_j.begin()+1, repr_j.end());
                return_amp = getValueByPath(edge.next_node->edges_.at(row_major_index), c_numb, repr_slice_i,repr_slice_j);

            }
            complex_number.returnToCache(c_numb);
            return return_amp;
        }



        CVec getVector(const vEdge& edge) {

            const std::size_t dim = std::accumulate(registers_sizes.begin(),registers_sizes.end(), 1, std::multiplies<>());
            // allocate resulting vector
            auto vec = CVec(dim, {0.0, 0.0});
            getVector(edge, Complex::one, 0, vec, dim);
            return vec;
        }


        void getVector(const vEdge& edge, const Complex& amp, std::size_t i, CVec& vec, std::size_t next) {
            // calculate new accumulated amplitude
            auto c_numb = complex_number.mulCached(edge.weight, amp);

            // base case
            if (edge.isTerminal()) {
                vec.at(i) = {CTEntry::val(c_numb.real), CTEntry::val(c_numb.img)};
                complex_number.returnToCache(c_numb);
                return;
            }

            auto offset = (next-i)/edge.next_node->edges_.size();

            for( auto k=0L; k < edge.next_node->edges_.size(); k++ ){
                if(!edge.next_node->edges_.at(k).weight.approximatelyZero()){
                    getVector(edge.next_node->edges_[k], c_numb, i+(k*offset), vec, i+((k+1)*offset));
                }
            }

            complex_number.returnToCache(c_numb);
        }



        static std::vector<unsigned long> getTreeNodes(const vEdge& edge){
            //MIXED BASIS DECODING
            std::vector<unsigned long > nodes_in_tree_from_edge;
            for (QuantumRegister j = edge.next_node->var_indx; j >= 0; j--) {
                nodes_in_tree_from_edge.push_back(j);
            }

            return nodes_in_tree_from_edge;
        }
        std::vector<unsigned long > getReprOfIndex(const std::vector<unsigned long> & nodes, const unsigned long i){
            std::vector<unsigned long> repr;
            repr.clear();
            // get representation
            auto quotient = i;
            auto remainder = 0UL;

            for (const unsigned long& index : nodes) {
                remainder = quotient % registers_sizes.at(index);
                quotient = quotient/registers_sizes.at(index);
                repr.push_back(remainder);
            }
            return repr;
        }


        void printVector(const vEdge& edge) {
            unsigned long long num_entries = std::accumulate(registers_sizes.begin(),registers_sizes.end(), 1, std::multiplies<>());

            auto nodes_in_tree = getTreeNodes(edge);

            for (auto i = 0ULL; i < num_entries; i++) {

                auto repr_i = getReprOfIndex(nodes_in_tree, i);
                //get amplitude
                const auto amplitude = getValueByPath(edge, repr_i);

                std::reverse(repr_i.begin(), repr_i.end());
                for(const unsigned long & coeff : repr_i ) {
                    std::cout << coeff;
                }
                repr_i.clear();

                constexpr auto precision = 3;
                // set fixed width to maximum of a printed number
                // (-) 0.precision plus/minus 0.precision i
                constexpr auto width = 1 + 2 + precision + 1 + 2 + precision + 1;
                std::cout << ": " << std::setw(width) << ComplexValue::toString(amplitude.r, amplitude.i, false, precision) << "\n";
            }
            std::cout << std::flush;
        }


    private:

        // check whether node represents a symmetric matrix or the identity
        void checkSpecialMatrices(mNode* node) {

            if (node->var_indx == -1)
                return;

            node->identity = false; // assume not identity
            node->symmetric  = false; // assume symmetric

            // check if matrix is symmetric
            auto number_of_edges = node->edges_.size();

            auto basic_dim = registers_sizes.at(node->var_indx);

            for(auto i = 0l; i< basic_dim;i++ ){
                if(!node->edges_.at(i*basic_dim + i).next_node->symmetric){
                    return;
                }
            }
            //TODO WHY RETURN IF DIAGONAL IS SYMMETRIC??
            //if (!node->edges_.at(0).next_node->symmetric || !node->edges_.at(3).next_node->symmetric) return;


            for( auto i = 0l; i<basic_dim; i++ ){
                for( auto j = 0l; j<basic_dim; j++ ){
                    if( i!=j ){
                        //row major indexing - enable optimization here
                        if( transpose(node->edges_.at(i*basic_dim + j)) != node->edges_.at(j*basic_dim + i) ){
                            return;
                        }
                    }
                }
            }
            //if (transpose(node->edges_.at(1)) != node->edges_.at(2)) return;

            node->symmetric = true;

            // check if matrix resembles identity
            for( auto i = 0l; i<basic_dim; i++ ){
                for( auto j = 0l; j<basic_dim; j++ ){
                        //row major indexing - enable optimization here
                        if(i==j){
                            if( !(node->edges_[i*basic_dim + j].next_node->identity) || (node->edges_[i*basic_dim + j].weight) != Complex::one ) return;
                        }
                        else{
                            if( (node->edges_[i*basic_dim + j].weight) != Complex::zero ) return;
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
            if (edge.next_node == nullptr || edge.isTerminal() || edge.next_node->symmetric) {
                return edge;
            }

            // check in compute table
            auto result = matrixTranspose.lookup(edge);
            if (result.next_node != nullptr) {
                return result;
            }

            std::vector<mEdge> new_edge{};
            auto basic_dim = registers_sizes.at(edge.next_node->var_indx);

            // transpose sub-matrices and rearrange as required
            for (auto i = 0U; i < basic_dim; i++) {
                for (auto j = 0U; j < basic_dim; j++) {
                    new_edge.at(basic_dim * i + j) = transpose(edge.next_node->edges_.at(basic_dim * j + i));
                }
            }
            // create new top node
            result = makeDDNode(edge.next_node->var_indx, new_edge);
            // adjust top weight
            auto c = complex_number.getTemporary();
            ComplexNumbers::mul(c, result.weight, edge.weight);
            result.weight = complex_number.lookup(c);

            // put in compute table
            matrixTranspose.insert(edge, result);
            return result;
        }
        mEdge conjugateTranspose(const mEdge& edge) {
            if (edge.next_node == nullptr)
                return edge;
            if (edge.isTerminal()) { // terminal case
                auto result = edge;
                result.weight    = ComplexNumbers::conj(edge.weight);
                return result;
            }

            // check if in compute table
            auto result = conjugateMatrixTranspose.lookup(edge);
            if (result.next_node != nullptr) {
                return result;
            }

            std::vector<mEdge> new_edge{};
            auto basic_dim = registers_sizes.at(edge.next_node->var_indx);

            // conjugate transpose submatrices and rearrange as required
            for (auto i = 0U; i < basic_dim; ++i) {
                for (auto j = 0U; j < basic_dim; ++j) {
                    new_edge.at(basic_dim * i + j)  = conjugateTranspose(edge.next_node->edges_.at(basic_dim * j + i));
                }
            }
            // create new top node
            result = makeDDNode(edge.next_node->var_indx, new_edge);

            auto c = complex_number.getTemporary();
            // adjust top weight including conjugate
            ComplexNumbers::mul(c, result.weight, ComplexNumbers::conj(edge.weight));
            result.weight = complex_number.lookup(c);


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