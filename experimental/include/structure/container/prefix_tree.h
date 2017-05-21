#pragma once

#include "algorithm/collections/collections.h"

#include <stdint.h>
#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * A prefix tree (Trie).
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 *   foreach: `O(n)`
 *   query: `O(l)`
 *   update: `O(l)`
 * Where `n` is the total number of characters and `l` is the length of the word.
 *
 * Note, word indices are 1-based.
 */
template<int ALPHABET_SIZE = 26, typename INDEX_T = int, typename ALPHA_T = uint8_t, typename WORD_T = std::string>
class prefix_tree {
public:
	typedef INDEX_T index_t;
	typedef ALPHA_T alpha_t;
    typedef WORD_T word_t;

protected:
	struct node_t {
        index_t word_id = 0;
        index_t num_ch = 0;
        alpha_t ord = 0;
        node_t* parent = nullptr;
        node_t *next[ALPHABET_SIZE];
        node_t() {
            for (int i = 0; i < ALPHABET_SIZE; i++) next[i] = nullptr;
        }
	};

    std::allocator<node_t> alloc;
    index_t nodes;
    node_t* root;
    std::vector<node_t*> words;

public:
	prefix_tree() : nodes(0), root(nullptr) {
        _init();
	}

    ~prefix_tree() {
        words.clear();
        _free_all(root);
    }

    void clear() {
        words.clear();
        _free_all(root);
        _init();
    }

    // Returns the 1-based index of the newly inserted word.
    template<typename It, typename F>
    index_t add(It begin, It end, F ordinal) {
        return _insert(root, begin, end, ordinal)->word_id;
    }
    template<typename F>
    index_t add(const word_t& word, F ordinal) {
        return add(word.begin(), word.end(), ordinal);
    }

    // Creates a new word by appending to an existing word.
    // The existing word remains unchanged in the tree.
    // Returns the 1-based index of the newly inserted word.
    template<typename It, typename F>
    index_t append(index_t word_id, It begin, It end, F ordinal) {
        return _insert(words[word_id], begin, end, ordinal)->word_id;
    }
    template<typename F>
    index_t append(index_t word_id, const word_t& word, F ordinal) {
        return append(word_id, word.begin(), word.end(), ordinal);
    }

    // Returns the 1-based index of the word if it exists, 0 otherwise.
    template<typename It, typename F>
    index_t find(It begin, It end, F ordinal) const {
        return _find(root, begin, end, ordinal)->word_id;
    }
    template<typename F>
    index_t find(const word_t& word, F ordinal) const {
        return find(word.begin(), word.end(), ordinal);
    }

    // Returns a copy of the word at the given 1-based index.
    template<typename F>
    word_t get(index_t word_id, F letter) const {
        vector<decltype(letter(0))> w;
        for (node_t* t = words[word_id]; t != root; t = t->parent) {
            w.push_back(letter(t->ord));
        }
        return word_t(w.rbegin(), w.rend());
    }

    // Iterates over all the words starting with the given prefix.
    template<typename V, typename F>
    void for_each(index_t word_id, V visitor, F letter) const {
        word_t w = get(word_id, letter);
        _for_each(words[word_id], w, visitor, letter);
    }
    template<typename V, typename F>
    void for_each(V visitor, F letter) const {
        for_each(0, visitor, letter);
    }

    // Removes the word at the given 1-based index.
    // Note: this causes the last word to be moved at the removed index.
    index_t erase(index_t word_id) {
        return _erase(words[word_id]);
    }
    template<typename It, typename F>
    index_t erase(It begin, It end, F ordinal) {
        return erase(find(begin, end, ordinal));
    }
    template<typename F>
    index_t erase(const word_t& word, F ordinal) {
        return erase(find(word, ordinal));
    }

    // Returns the number of words.
    index_t num_words() const {
        return (index_t)words.size() - 1;
    }

    // Returns the number of letters in the tree.
    // Letters in shared prefixes are counted only once.
    index_t num_letters() const {
        return nodes - 1;
    }

    // helper functors for ordinal <-> letter conversion
	static alpha_t ordinal_digit(char c) { return alpha_t(c - '0'); }
    static alpha_t ordinal_lower_alpha(char c) { return alpha_t(c - 'a'); }
	static alpha_t ordinal_upper_alpha(char c) { return alpha_t(c - 'A'); }
    static char letter_digit(alpha_t o) { return char(o + '0'); }
    static char letter_lower_alpha(alpha_t o) { return char(o + 'a'); }
    static char letter_upper_alpha(alpha_t o) { return char(o + 'A'); }

private:
    template<typename It, typename F>
    node_t* _find(node_t* t, It begin, It end, F ordinal) const {
        for (auto it = begin; it != end; ++it) {
            alpha_t o = ordinal(*it);
            if (!t->next[o]) return root;
            t = t->next[o];
        }
        return t;
    }

    template<typename It, typename F>
    node_t* _insert(node_t* t, It begin, It end, F ordinal) {
        for (auto it = begin; it != end; ++it) {
            alpha_t o = ordinal(*it);
            if (!t->next[o]) {
                t->next[o] = _buy();
                t->next[o]->parent = t;
                t->next[o]->ord = o;
                t->num_ch++;
            }
            t = t->next[o];
        }
        if (t != root && !t->word_id) {
            t->word_id = (index_t)words.size();
            words.push_back(t);
        }
        return t;
    }

    index_t _erase(node_t* t) {
        index_t word_id = t->word_id;
        if (word_id == 0) return 0;
        swap(words[word_id], words.back());
        words[word_id]->word_id = word_id;
        words.pop_back();
        t->word_id = 0;
        while (t != root && t->num_ch == 0 && t->word_id == 0) {
            node_t* p = t->parent;
            p->num_ch--;
            p->next[t->ord] = nullptr;
            _free(t);
            t = p;
        }
        return word_id;
    }

    template<typename V, typename F>
    void _for_each(node_t* t, word_t& w, V visitor, F letter) const {
        if (t->word_id != 0) visitor(w, t->word_id);
        for (int o = 0; o < ALPHABET_SIZE; o++) {
            if (t->next[o] == nullptr) continue;
            w.push_back(letter(o));
            _for_each(t->next[o], w, visitor, letter);
            w.pop_back();
        }
    }

    void _init() {
        root = _buy();
        words.push_back(root);
    }

    node_t* _buy() {
        node_t* t = alloc.allocate(1);
        alloc.construct(t);
        nodes++;
        return t;
    }

    void _free(node_t* t) {
        nodes--;
        alloc.destroy(t);
        alloc.deallocate(t, 1);
    }

    void _free_all(node_t* t) {
        if (t == nullptr) return;
        for (int o = 0; o < ALPHABET_SIZE; o++) {
            _free_all(t->next[o]);
        }
        _free(t);
    }
};

} // container
} // altruct
