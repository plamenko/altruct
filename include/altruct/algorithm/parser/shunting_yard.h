#pragma once

#include <vector>
#include <functional>

namespace altruct {
namespace parser {
namespace shunting_yard {

/**
 * Infix token.
 */
struct infix_token {
    enum class type_t { OPERAND, OPERATOR, LEFT_PARENTHESIS, RIGHT_PARENTHESIS, FUNCTION, SEPARATOR };

    type_t type;
    int id;
};

/**
 * Postfix token.
 */
struct postfix_token {
    enum class type_t { OPERAND, FUNCTION };

    type_t type;
    int id;
    int num_args;

    bool operator == (const postfix_token& rhs) const {
        return type == rhs.type && id == rhs.id && num_args == rhs.num_args;
    }
};

/**
 * Infix operator descriptor.
 *
 * Operator with the bigger precedence value has priority.
 * In case the precedence values are the same, associativity determines priority.
 */
struct operator_desc {
    enum class assoc_t { LEFT, RIGHT };
    enum class arity_t { UNARY, BINARY };

    int precedence;
    assoc_t associativity;
    arity_t arity;

    bool operator < (const operator_desc &rhs) const {
        if (precedence != rhs.precedence) return (precedence < rhs.precedence);
        return (associativity == assoc_t::LEFT);
    }

    int num_args() const {
        return (arity == arity_t::UNARY) ? 1 : 2;
    }
};

/**
 * Converts expresion from infix to postfix (a.k.a. Reverse-Polish) notation.
 *
 * Supports:
 *   + unary and binary operators
 *   + operator precedence and asociativity
 *   + parentheses
 *   + functions with any number of arguments
 *   + multiple expressions (separated by comma)
 *
 * Complexity: O(n)
 *
 * param infix_tokens - infix tokens to convert
 * param operators    - operator descriptors
 */
template<typename T = void>
std::vector<postfix_token> infix_to_postfix(
      std::vector<infix_token> &infix_tokens,
      std::vector<operator_desc> operators) {
    int num_args = 0;
    std::vector<postfix_token> res;
    std::vector<infix_token> stk;
    auto pop_operators = [&](int id) {
        while (!stk.empty() && stk.back().type == infix_token::type_t::OPERATOR &&
            (id < 0 || operators[id] < operators[stk.back().id])) {
            num_args = operators[stk.back().id].num_args();
            res.push_back({ postfix_token::type_t::FUNCTION, stk.back().id, num_args });
            stk.pop_back();
        }
    };
    for (const auto& t : infix_tokens) {
        switch (t.type) {
            case infix_token::type_t::OPERAND:
                res.push_back({ postfix_token::type_t::OPERAND, t.id, 0 });
                break;
            case infix_token::type_t::OPERATOR:
                pop_operators(t.id);
                stk.push_back(t);
                break;
            case infix_token::type_t::LEFT_PARENTHESIS:
                stk.push_back(t);
                break;
            case infix_token::type_t::RIGHT_PARENTHESIS:
                pop_operators(-1);
                num_args = 1;
                while (!stk.empty() && stk.back().type == infix_token::type_t::SEPARATOR) {
                    num_args++; stk.pop_back();
                }
                if (!stk.empty() && stk.back().type == infix_token::type_t::LEFT_PARENTHESIS) {
                    stk.pop_back();
                }
                if (!stk.empty() && stk.back().type == infix_token::type_t::FUNCTION) {
                    res.push_back({ postfix_token::type_t::FUNCTION, stk.back().id, num_args });
                    stk.pop_back();
                }
                break;
            case infix_token::type_t::FUNCTION:
                stk.push_back(t);
                break;
            case infix_token::type_t::SEPARATOR:
                pop_operators(-1);
                stk.push_back(t);
                break;
        }
    }
    pop_operators(-1);
    return res;
}

/**
 * Evaluates expression in postfix (a.k.a. Reverse-Polish) notation.
 *
 * Complexity: O(n)
 *
 * param postfix_tokens  - expression to evaluate given as postfix tokens
 * param operand_valuess - values of the operands used in the given expression
 * param evaluators      - evaluators used to evaluate the functions used in the given expression
 */
template<typename T, typename VE>
std::vector<T> evaluate_postfix(
      const std::vector<postfix_token> &postfix_tokens,
      const std::vector<T> &operand_values,
      const VE& evaluators) {
    std::vector<T> stk;
    for (const auto& t : postfix_tokens) {
        switch (t.type) {
            case postfix_token::type_t::OPERAND:
                stk.push_back(operand_values[t.id]);
                break;
            case postfix_token::type_t::FUNCTION:
                int pos = (int)stk.size() - t.num_args;
                T r = evaluators[t.id](&stk[pos]);
                stk.resize(pos);
                stk.push_back(r);
                break;
        }
    }
    return stk;
}

/**
 * Example implementation for basic math operators.
 */
template<typename T>
struct basic_math {
    enum evaluators { NEG, MUL, DIV, ADD, SUB };

    static const std::vector<operator_desc>& operators() {
        static std::vector<operator_desc> operators{
            { 20, operator_desc::assoc_t::RIGHT, operator_desc::arity_t::UNARY },  // NEG
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // MUL
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // DIV
            { 10, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // ADD
            { 10, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // SUB
        };
        return operators;
    }

    static const std::vector<std::function<T(T*)>>& evaluators() {
        static std::vector<std::function<T(T*)>> evaluators{
            [](T* a){ return -a[0]; },       // NEG
            [](T* a){ return a[0] * a[1]; }, // MUL
            [](T* a){ return a[0] / a[1]; }, // DIV
            [](T* a){ return a[0] + a[1]; }, // ADD
            [](T* a){ return a[0] - a[1]; }, // SUB
        };
        return evaluators;
    }
};

/**
 * Example implementation for integer math operators and functions.
 */
template<typename T>
struct integer_math {
    enum evaluators {
        NEG, MUL, DIV, MOD, ADD, SUB,
        SHL, SHR, ROL, ROR,
        LTE, GTE, LT, GT, EQ, NEQ,
        BIT_NOT, BIT_AND, BIT_XOR, BIT_OR,
        LOG_NOT, LOG_AND, LOG_XOR, LOG_OR,
        SQRT, POW,
    };

    static const std::vector<operator_desc>& operators() {
        static std::vector<operator_desc> operators{
            // arithmetic
            { 20, operator_desc::assoc_t::RIGHT, operator_desc::arity_t::UNARY },  // NEG
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // MUL
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // DIV
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // MOD
            { 10, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // ADD
            { 10, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // SUB
            // shifts
            { 9, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // SHL
            { 9, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // SHR
            { 9, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // ROL
            { 9, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // ROR
            // relational
            { 8, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // LTE
            { 8, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // GTE
            { 8, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // LT
            { 8, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // GT
            { 7, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // EQ
            { 7, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // NEQ
            // bitwise
            { 20, operator_desc::assoc_t::RIGHT, operator_desc::arity_t::UNARY },  // BIT_NOT
            { 6, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // BIT_AND
            { 5, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // BIT_XOR
            { 4, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // BIT_OR
            // logical
            { 20, operator_desc::assoc_t::RIGHT, operator_desc::arity_t::UNARY },  // LOG_NOT
            { 3, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },   // LOG_AND
            { 2, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },   // LOG_XOR
            { 1, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },   // LOG_OR
        };
        return operators;
    }

    static const std::vector<std::function<T(T*)>>& evaluators() {
        static std::vector<std::function<T(T*)>> evaluators{
            // arithmetic
            [](T* a){ return -a[0]; },       // NEG
            [](T* a){ return a[0] * a[1]; }, // MUL
            [](T* a){ return a[0] / a[1]; }, // DIV
            [](T* a){ return a[0] % a[1]; }, // MOD
            [](T* a){ return a[0] + a[1]; }, // ADD
            [](T* a){ return a[0] - a[1]; }, // SUB
            // shifts
            [](T* a){ return a[0] << a[1]; }, // SHL
            [](T* a){ return a[0] >> a[1]; }, // SHR
            [](T* a){ return a[0] << a[1]; }, // ROL // TODO
            [](T* a){ return a[0] >> a[1]; }, // ROR // TODO
            // relational
            [](T* a){ return a[0] <= a[1]; }, // LTE
            [](T* a){ return a[0] >= a[1]; }, // GTE
            [](T* a){ return a[0] < a[1]; },  // LT
            [](T* a){ return a[0] > a[1]; },  // GT
            [](T* a){ return a[0] == a[1]; }, // EQ
            [](T* a){ return a[0] != a[1]; }, // NEQ
            // bitwise
            [](T* a){ return ~a[0]; },        // BIT_NOT
            [](T* a){ return a[0] & a[1]; },  // BIT_AND
            [](T* a){ return a[0] ^ a[1]; },  // BIT_XOR
            [](T* a){ return a[0] | a[1]; },  // BIT_OR
            // logical
            [](T* a){ return !a[0]; },         // LOG_NOT
            [](T* a){ return a[0] && a[1]; },  // LOG_AND
            [](T* a){ return !!a[0] ^ !!a[1]; },  // LOG_XOR
            [](T* a){ return a[0] || a[1]; },  // LOG_OR
            // functions
            [](T* a){ return sqrtT(a[0]); },      // SQRT
            [](T* a){ return powT(a[0], a[1]); }, // POW
        };
        return evaluators;
    }
};

/**
 * Example implementation for floating-point math operators and functions.
 */
template<typename T>
struct floating_point_math {
    enum evaluators {
        NEG, MUL, DIV, ADD, SUB,
        SQRT, POW, EXP, LOG,
        SIN, COS, TAN,
        SINH, COSH, TANH,
        ASIN, ACOS, ATAN, ATAN2,
        ASINH, ACOSH, ATANH
    };

    static const std::vector<operator_desc>& operators() {
        static std::vector<operator_desc> operators{
            // arithmetic
            { 20, operator_desc::assoc_t::RIGHT, operator_desc::arity_t::UNARY },  // NEG
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // MUL
            { 11, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // DIV
            { 10, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // ADD
            { 10, operator_desc::assoc_t::LEFT, operator_desc::arity_t::BINARY },  // SUB
        };
        return operators;
    }

    static const std::vector<std::function<T(T*)>>& evaluators() {
        static std::vector<std::function<T(T*)>> evaluators{
            [](T* a){ return -a[0]; },       // NEG
            [](T* a){ return a[0] * a[1]; }, // MUL
            [](T* a){ return a[0] / a[1]; }, // DIV
            [](T* a){ return a[0] + a[1]; }, // ADD
            [](T* a){ return a[0] - a[1]; }, // SUB

            [](T* a){ return sqrt(a[0]); },        // SQRT
            [](T* a){ return pow(a[0], a[1]); },   // POW
            [](T* a){ return exp(a[0]); },         // EXP
            [](T* a){ return log(a[0]); },         // LOG

            [](T* a){ return sin(a[0]); },         // SIN
            [](T* a){ return cos(a[0]); },         // COS
            [](T* a){ return tan(a[0]); },         // TAN

            [](T* a){ return sinh(a[0]); },        // SINH
            [](T* a){ return cosh(a[0]); },        // COSH
            [](T* a){ return tanh(a[0]); },        // TANH

            [](T* a){ return asin(a[0]); },        // ASIN
            [](T* a){ return acos(a[0]); },        // ACOS
            [](T* a){ return atan(a[0]); },        // ATAN
            [](T* a){ return atan2(a[0], a[1]); }, // ATAN2

            [](T* a){ return asinh(a[0]); },       // ASINH
            [](T* a){ return acosh(a[0]); },       // ACOSH
            [](T* a){ return atanh(a[0]); },       // ATANH
        };
        return evaluators;
    }
};

} // shunting_yard
} // parser
} // altruct
