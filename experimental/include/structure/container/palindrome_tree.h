#pragma once

#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * A tree that enumerates palindrome substrings.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 */
template<int ALPHABET_SIZE = 26>
class palindrome_tree {
public:
	typedef char ordinal_t;

	// node that represents a palindromic substring
	struct node_t {
		size_t len;      // length of this palindromic substring
		size_t pos;      // position of the first occurence of this palindromic substring within the string
		size_t cnt;      // multiplicity of this palindromic substring
		size_t depth;    // depth in the suffix chain of this node
		size_t suff;     // node-index of the largest palindromic suffix of this node
		size_t next[ALPHABET_SIZE]; // "A".next['x'] --> "xAx"
	};

	static const size_t NIL = 0;
	static const size_t NEGAT = 1;
	static const size_t EMPTY = 2;
	static const size_t RESERVED = 3;

	std::vector<ordinal_t> _string; // string of letter ordinals (e.g. 'a' is 0)
	std::vector<node_t> _nodes;     // node container
	size_t _suff;                   // node-index of the current longest palindromic suffix
	size_t _total;                  // total number of palindromic substrings (counting multiplicities)

	palindrome_tree() {
		_init();
	}

	template<typename It, typename F>
	palindrome_tree(It begin, It end, F ordinal) {
		_init();
		add_all(begin, end, ordinal);
	}

	void _init() {
		_string.clear();
		_nodes.resize(RESERVED);
		_nodes[NEGAT].len = -1; _nodes[NEGAT].suff = NEGAT;
		_nodes[EMPTY].len = 0; _nodes[EMPTY].suff = NEGAT;
		_suff = EMPTY;
		_total = 0;
	}

	template<typename It, typename F>
	size_t add_all(It begin, It end, F ordinal) {
		size_t c = 0;
		for (auto it = begin; it != end; ++it) {
			c += add(*it, ordinal);
		}
		return c;
	}

	template<typename T, typename F>
	size_t add(const T& t, F ordinal) {
		return add(ordinal(t));
	}
	
	size_t add(ordinal_t ord) {
		_string.push_back(ord);
		auto i = _find_suffix(_suff, ord);
		_suff = _nodes[i].next[ord];
		if (_suff != NIL) {
			_nodes[_suff].cnt++;
			_total += _nodes[_suff].depth;
			return 0;
		}
		_suff = _nodes.size();
		_nodes.push_back({});
		_nodes[i].next[ord] = _suff;
		auto suff2 = _find_suffix2(i, ord);
		_nodes[_suff].len = _nodes[i].len + 2;
		_nodes[_suff].pos = _string.size() - _nodes[_suff].len;
		_nodes[_suff].cnt = 1;
		_nodes[_suff].suff = suff2;
		_nodes[_suff].depth = _nodes[suff2].depth + 1;
		_total += _nodes[_suff].depth;
		return 1;
	}

	size_t _find_suffix2(size_t i, ordinal_t ord) {
		if (i == NEGAT) return EMPTY;
		i = _find_suffix(_nodes[i].suff, ord);
		return _nodes[i].next[ord];
	}

	size_t _find_suffix(size_t i, ordinal_t ord) {
		auto sz = _string.size();
		while (sz < _nodes[i].len + 2 || _string[sz - _nodes[i].len - 2] != ord) {
			i = _nodes[i].suff;
		}
		return i;
	}

	// This should be called only once after all elements are added!
	void propagate() {
		for (auto i = _nodes.size() - 1; i >= RESERVED; i--) {
			auto suff = _nodes[i].suff;
			_nodes[suff].cnt += _nodes[i].cnt;
		}
	}

	// Returns the number of total palindromic substrings, counting their multiplicities.
	size_t total() const {
		return _total;
	}

	// Returns the number of distinct palindromic substrings, each counted only once.
	size_t distinct() const {
		return _nodes.size() - RESERVED;
	}

	// Returns the index of the first node.
	size_t first() const {
		return RESERVED;
	}

	// Returns the number of nodes.
	size_t size() const {
		return _nodes.size();
	}

	// Accesses node by its index.
	node_t& operator[] (size_t index) {
		return _nodes[index];
	}
};

} // container
} // altruct
