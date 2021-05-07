#include <memory>
#include "token.h"

class Expr {
public:
    int ssa;
    Expr() {}
    virtual ~Expr() = default;
};

class BinOp : public Expr {
public:
    std::unique_ptr<Expr> left;
    std::unique_ptr<Operator> op;
    std::unique_ptr<Expr> right;
    BinOp(std::unique_ptr<Expr>&& _left, std::unique_ptr<Operator>&& _op, std::unique_ptr<Expr>&& _right):
        Expr(), left(std::move(_left)), op(std::move(_op)), right(std::move(_right)) {}
};

class UnaryOp : public Expr {
public:
    std::unique_ptr<Operator> op;
    std::unique_ptr<Expr> operand;
    UnaryOp(std::unique_ptr<Operator>&& _op, std::unique_ptr<Expr>&& _operand):
        Expr(), op(std::move(_op)), operand(std::move(_operand)) {}
};

class Literal : public Expr {
//private:
public:
    std::unique_ptr<Token> value;
public:
    Literal(std::unique_ptr<Token>&& _value): Expr(), value(std::move(_value)) {}
};
