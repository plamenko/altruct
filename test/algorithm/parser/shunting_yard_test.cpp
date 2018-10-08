#include "altruct/algorithm/parser/shunting_yard.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::parser;

TEST(shunting_yard_test, convert_and_evaluate_for_mod) {
    typedef modulo<int, 1000000007, modulo_storage::CONSTANT> mod;
    using tok = shunting_yard::infix_token::type_t;
    using math = shunting_yard::basic_math<mod>;
    vector<shunting_yard::infix_token> expr;
    expr.push_back({ tok::OPERAND, 0 });
    expr.push_back({ tok::OPERATOR, math::SUB });
    expr.push_back({ tok::OPERAND, 1 });
    expr.push_back({ tok::OPERATOR, math::DIV });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::OPERAND, 2 });
    expr.push_back({ tok::OPERATOR, math::ADD });
    expr.push_back({ tok::OPERAND, 3 });
    expr.push_back({ tok::OPERATOR, math::MUL });
    expr.push_back({ tok::OPERATOR, math::NEG });
    expr.push_back({ tok::OPERAND, 4 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::SEPARATOR });
    expr.push_back({ tok::OPERAND, 5 });
    expr.push_back({ tok::SEPARATOR });
    expr.push_back({ tok::OPERAND, 6 });
    expr.push_back({ tok::OPERATOR, math::SUB });
    expr.push_back({ tok::OPERATOR, math::NEG });
    expr.push_back({ tok::OPERAND, 7 });

    using ptok = shunting_yard::postfix_token::type_t;
    auto expected = vector<shunting_yard::postfix_token>{
        { ptok::OPERAND, 0, 0 },
        { ptok::OPERAND, 1, 0 },
        { ptok::OPERAND, 2, 0 },
        { ptok::OPERAND, 3, 0 },
        { ptok::OPERAND, 4, 0 },
        { ptok::FUNCTION, math::NEG, 1 },
        { ptok::FUNCTION, math::MUL, 2 },
        { ptok::FUNCTION, math::ADD, 2 },
        { ptok::FUNCTION, math::DIV, 2 },
        { ptok::FUNCTION, math::SUB, 2 },
        { ptok::OPERAND, 5, 0 },
        { ptok::OPERAND, 6, 0 },
        { ptok::OPERAND, 7, 0 },
        { ptok::FUNCTION, math::NEG, 1 },
        { ptok::FUNCTION, math::SUB, 2 },
    };
    auto rpn_expr = shunting_yard::infix_to_postfix<>(expr, math::operators());
    EXPECT_EQ(expected, rpn_expr);

    vector<mod> operands;
    operands.push_back(8);
    operands.push_back(7);
    operands.push_back(2);
    operands.push_back(3);
    operands.push_back(5);
    operands.push_back(-4);
    operands.push_back(6);
    operands.push_back(9);

    auto res = shunting_yard::evaluate_postfix(rpn_expr, operands, math::evaluators());
    EXPECT_EQ((vector<mod>{76923086, -4, 15}), res);
}

TEST(shunting_yard_test, convert_and_evaluate_for_int) {
    using tok = shunting_yard::infix_token::type_t;
    using math = shunting_yard::integer_math<int>;
    vector<shunting_yard::infix_token> expr;
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::OPERATOR, math::NEG });
    expr.push_back({ tok::FUNCTION, math::SQRT });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::FUNCTION, math::POW });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::OPERAND, 0 });
    expr.push_back({ tok::OPERATOR, math::ADD });
    expr.push_back({ tok::OPERATOR, math::NEG });
    expr.push_back({ tok::OPERAND, 1 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::OPERATOR, math::MUL });
    expr.push_back({ tok::OPERAND, 2 });
    expr.push_back({ tok::SEPARATOR });
    expr.push_back({ tok::OPERAND, 3 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::SEPARATOR });
    expr.push_back({ tok::OPERAND, 1 });
    expr.push_back({ tok::OPERATOR, math::DIV });
    expr.push_back({ tok::OPERAND, 0 });
    expr.push_back({ tok::OPERATOR, math::SUB });
    expr.push_back({ tok::OPERAND, 3 });
    expr.push_back({ tok::OPERATOR, math::MOD });
    expr.push_back({ tok::OPERAND, 2 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    // TODO: shifts, relational, bitwise, logical

    auto rpn_expr = shunting_yard::infix_to_postfix<>(expr, math::operators());
    // TODO: EXPECT_EQ

    vector<int> operands;
    operands.push_back(1);
    operands.push_back(2);
    operands.push_back(3);
    operands.push_back(4);

    auto res = shunting_yard::evaluate_postfix(rpn_expr, operands, math::evaluators());
    EXPECT_EQ((vector<int>{-9, 1}), res);
}

TEST(shunting_yard_test, convert_and_evaluate_for_double) {
    using tok = shunting_yard::infix_token::type_t;
    using math = shunting_yard::floating_point_math<double>;
    vector<shunting_yard::infix_token> expr;
    expr.push_back({ tok::FUNCTION, math::POW });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::FUNCTION, math::LOG });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::FUNCTION, math::EXP });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::OPERAND, 0 });
    expr.push_back({ tok::OPERATOR, math::ADD });
    expr.push_back({ tok::OPERATOR, math::NEG });
    expr.push_back({ tok::OPERAND, 1 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::OPERATOR, math::DIV });
    expr.push_back({ tok::OPERAND, 2 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::SEPARATOR });
    expr.push_back({ tok::OPERAND, 3 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    expr.push_back({ tok::SEPARATOR });
    expr.push_back({ tok::FUNCTION, math::SQRT });
    expr.push_back({ tok::LEFT_PARENTHESIS });
    expr.push_back({ tok::OPERAND, 3 });
    expr.push_back({ tok::OPERATOR, math::MUL });
    expr.push_back({ tok::OPERAND, 2 });
    expr.push_back({ tok::OPERATOR, math::SUB });
    expr.push_back({ tok::OPERAND, 0 });
    expr.push_back({ tok::RIGHT_PARENTHESIS });
    // TODO: trigonometric

    auto rpn_expr = shunting_yard::infix_to_postfix<>(expr, math::operators());
    // TODO: EXPECT_EQ

    vector<double> operands;
    operands.push_back(1);
    operands.push_back(2);
    operands.push_back(3);
    operands.push_back(4);

    auto res = shunting_yard::evaluate_postfix(rpn_expr, operands, math::evaluators());
    EXPECT_NEAR(1.0 / 81, res[0], 1e-9);
    EXPECT_NEAR(sqrt(11.0), res[1], 1e-9);
}
