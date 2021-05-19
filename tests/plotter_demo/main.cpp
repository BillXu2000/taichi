#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <stack>
#include <set>
#include <cassert>
#include "expr.h"

#include "taichi/ir/ir_builder.h"
#include "taichi/ir/statements.h"
#include "taichi/program/program.h"
#include "taichi/ir/ir.h"

#include "taichi/gui/gui.h"

std::vector<std::unique_ptr<Token> > get_lexical(std::string plain, std::set<std::string> &names, std::set<std::string> reserved_funcs) {
    using namespace std;
    vector<unique_ptr<Token> > lex;
    auto is_number_or_dot = [](char ch) {
        return isdigit(ch) || ch == '.';
    };
    bool expect_operator = false;
    lex.push_back(make_unique<OpenParen>());
    for (int i = 0; plain[i];) {
        if (plain[i] == ';') break;
        if (isspace(plain[i])) i++;
        else if (is_number_or_dot(plain[i])) {
            assert(!expect_operator);
            expect_operator = true;
            int j;
            for (j = i + 1; is_number_or_dot(plain[j]); j++);
            lex.push_back(make_unique<Constant>(atof(plain.substr(i, j - i).c_str())));
            i = j;
        }
        else if (isalpha(plain[i]) || plain[i] == '_') {
            assert(!expect_operator);
            int j;
            for (j = i + 1; isalpha(plain[j]) || plain[j] == '_'; j++);
            string name = plain.substr(i, j - i);
            if (reserved_funcs.find(name) == reserved_funcs.end()) {
                lex.push_back(make_unique<Variable>(name));
                names.insert(name);
                expect_operator = true;
            }
            else lex.push_back(make_unique<UnaryToken>(name));
            i = j;
        }
        else {
            if (plain[i] == '-') {
                if (expect_operator) lex.push_back(make_unique<Sub>());
                else lex.push_back(make_unique<UnaryToken>("-"));
                expect_operator = false;
            }
            else if (plain[i] == '(') {
                assert(!expect_operator);
                lex.push_back(make_unique<OpenParen>());
            }
            else if (plain[i] == ')') {
                assert(expect_operator);
                lex.push_back(make_unique<CloseParen>());
            }
            else {
                assert(expect_operator);
                expect_operator = false;
                switch (plain[i]) {
                    case '+':
                        lex.push_back(make_unique<Add>());
                        break;
                    case '*':
                        lex.push_back(make_unique<Mult>());
                        break;
                    case '/':
                        lex.push_back(make_unique<Div>());
                        break;
                }
            }
            i++;
        }
    }
    lex.push_back(make_unique<CloseParen>());
    return lex;
}

Expr* get_ast(std::vector<std::unique_ptr<Token> >& lex) {
    using namespace std;
    stack<unique_ptr<Operator> > ops;
    stack<unique_ptr<Expr> > nodes;
    for (auto &i : lex) {
        auto op = dynamic_cast<Operator*>(i.get());
        if (op) {
            if (dynamic_cast<OpenParen*>(op) || dynamic_cast<UnaryToken*>(op)) {
                i.release();
                ops.push(unique_ptr<Operator>(op));
                continue;
            }
            for (;!ops.empty() && ops.top()->get_precedence() >= op->get_precedence();) {
                if (dynamic_cast<OpenParen*>(ops.top().get())) {
                    ops.pop();
                    break;
                }
                auto right = nodes.top().get();
                nodes.top().release();
                nodes.pop();
                if (dynamic_cast<UnaryToken*>(ops.top().get())) {
                    nodes.push(make_unique<UnaryOp>(move(ops.top()), unique_ptr<Expr>(right)));
                }
                else {
                    auto left = nodes.top().get();
                    nodes.top().release();
                    nodes.pop();
                    nodes.push(make_unique<BinOp>(unique_ptr<Expr>(left), move(ops.top()), unique_ptr<Expr>(right)));
                }
                ops.pop();
            }
            if (!dynamic_cast<CloseParen*>(op)) {
                i.release();
                ops.push(unique_ptr<Operator>(op));
            }
        }
        else nodes.push(make_unique<Literal>(move(i)));
    }
    auto ans = nodes.top().get();
    nodes.top().release();
    return ans;
}

void print(std::unique_ptr<Expr> &i) {
    using namespace std;
    auto p = i.get();
    auto bin = dynamic_cast<BinOp*>(p);
    auto unary = dynamic_cast<UnaryOp*>(p);
    if (bin) {
        cout << bin->op->to_string() << "(";
        print(bin->left);
        cout << ", ";
        print(bin->right);
        cout << ")";
    }
    else if (unary) {
        cout << unary->op->to_string() << "(";
        print(unary->operand);
        cout << ")";
    }
    else cout << dynamic_cast<Literal*>(p)->value->to_string();
}

void init_IR(std::fstream &fs, std::set<std::string> &names) {
    fs << ""
"import taichi as ti\n"
//"n_IR = 100000\n"
"real = ti.f32\n";
//"ssa = ti.field(dtype = real, shape = n_IR, needs_grad = True)\n";
//"ans = ti.field(dtype = real, shape = (), needs_grad = True)\n";
    for (auto i : names) {
        if (i != "x" && i != "y") fs << i + " = ti.field(dtype = real, shape = (), needs_grad = True)\n";
    }
    fs << ""
"@ti.func\n"
"def generated_ir(x, y):\n";
}

void dummy_IR(std::unique_ptr<Expr> &i, std::fstream &fs, int &cnt) {
    using namespace std;
    auto p = i.get();
    auto bin = dynamic_cast<BinOp*>(p);
    auto unary = dynamic_cast<UnaryOp*>(p);
    if (bin) {
        dummy_IR(bin->left, fs, cnt);
        dummy_IR(bin->right, fs, cnt);
        p->ssa = cnt++;
        //fs << "\tssa[" + to_string(bin->ssa) + "] = ssa[" + to_string(bin->left->ssa) + "] " + bin->op->to_string() + " ssa[" + to_string(bin->right->ssa) + "]\n";
        fs << "\tssa" + to_string(bin->ssa) + " = ssa" + to_string(bin->left->ssa) + " " + bin->op->to_string() + " ssa" + to_string(bin->right->ssa) + "\n";
    }
    else if (unary) {
        dummy_IR(unary->operand, fs, cnt);
        p->ssa = cnt++;
        //fs << "\tssa[" + to_string(unary->ssa) + "] = " + unary->op->to_string() + "(ssa[" + to_string(unary->operand->ssa) + "])\n";
        fs << "\tssa" + to_string(unary->ssa) + " = " + unary->op->to_string() + "(ssa" + to_string(unary->operand->ssa) + ")\n";
    }
    else {
        p->ssa = cnt++;
        //fs << "\tssa[" + to_string(p->ssa) + "] = " + dynamic_cast<Literal*>(p)->value->to_string() + "\n";
        fs << "\tssa" + to_string(p->ssa) + " = " + dynamic_cast<Literal*>(p)->value->to_string() + "\n";
    }
}

taichi::lang::Stmt* build_IR(std::unique_ptr<Expr> &i, taichi::lang::IRBuilder &builder, std::map<std::string, taichi::lang::Stmt*> para) {
    using namespace taichi::lang;
    auto p = i.get();
    auto bin = dynamic_cast<BinOp*>(p);
    auto unary = dynamic_cast<UnaryOp*>(p);
    if (bin) {
        auto *lhs = build_IR(bin->left, builder, para);
        auto *rhs = build_IR(bin->right, builder, para);
        auto *op = bin->op.get();
        if (dynamic_cast<Add*>(op)) return builder.create_add(lhs, rhs);
        if (dynamic_cast<Sub*>(op)) return builder.create_sub(lhs, rhs);
        if (dynamic_cast<Mult*>(op)) return builder.create_mul(lhs, rhs);
        if (dynamic_cast<Div*>(op)) return builder.create_div(lhs, rhs);
    }
    else if (unary) {
        auto *operand = build_IR(unary->operand, builder, para);
        auto *op = dynamic_cast<UnaryToken*>(unary->op.get());
        if (op->name == "-") return builder.create_neg(operand);
        if (op->name == "sin") return builder.create_sin(operand);
        if (op->name == "cos") return builder.create_cos(operand);
    }
    else {
        std::string name = dynamic_cast<Literal*>(p)->value->to_string();
        if (isdigit(name[0]) || name[0] == '-') return builder.get_float32(atof(name.c_str()));
        return para[name];
    }
    return nullptr;
}

void run_builder(taichi::lang::IRBuilder &builder) {
    using namespace taichi::lang;
    auto prog_ = Program(arch_from_name("x64"));
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    prog_.config = config_print_ir;
    prog_.materialize_layout();
    auto ker = std::make_unique<Kernel>(prog_, builder.extract_ir());
    auto launch_ctx = ker->make_launch_context();
    /*ker->insert_arg(PrimitiveType::f32, false);
    ker->insert_arg(PrimitiveType::f32, false);
    std::fstream fs("./512.out", std::fstream::out);
    const int N = 512;
    for (int i = 0; i < N; i++) {
        double x = double(i) / N;
        for (int j = 0; j < N; j++) {
            double y = double(j) / N;
            launch_ctx.set_arg_float(0, x);
            launch_ctx.set_arg_float(1, y);
            (*ker)(launch_ctx);
            fs << prog_.fetch_result<float>(0) << "\n";
        }
    }*/
    /*launch_ctx.w_float(0, 2);
    launch_ctx.set_arg_float(1, 3);
    (*ker)(launch_ctx);
    fs << prog_.fetch_result<float>(0) << "\n";*/
    (*ker)(launch_ctx);
}

namespace taichi {
namespace lang {
void test_snode() {
    auto prog_ = Program(arch_from_name("x64"));
    //prog_.thread_pool = std::make_unique<ThreadPool>(1);
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    prog_.config = config_print_ir;  // print ir
    auto root = prog_.snode_root.get();
    auto &dense = root->dense(Index(0), 10);
    auto &place = dense.insert_children(SNodeType::place);
    place.dt = PrimitiveType::i32;
    prog_.materialize_layout();


    IRBuilder builder;
    auto *zero = builder.get_int32(0);
    auto *ten = builder.get_int32(10);
    auto *loop = builder.create_range_for(zero, ten, 1, 0, 4);
    {
        builder.set_insertion_point_to_loop_begin(loop);
        auto *index = builder.get_loop_index(loop, 0);
        auto *ptr = builder.insert(Stmt::make<GlobalPtrStmt>(LaneAttribute<SNode *>(&place), std::vector<Stmt*>(1, index)));
        builder.insert(std::make_unique<GlobalStoreStmt>(ptr, index));
        builder.set_insertion_point_to_after(loop);
    }
    auto *index = builder.get_int32(4);
    auto *ret_ptr = builder.insert(Stmt::make<GlobalPtrStmt>(LaneAttribute<SNode *>(&place), std::vector<Stmt*>(1, index)));
    auto *ret_val = builder.insert(std::make_unique<GlobalLoadStmt>(ret_ptr));
    builder.create_return(ret_val);


    auto block = builder.extract_ir();
    auto ker = std::make_unique<Kernel>(prog_, std::move(block));
    auto launch_ctx = ker->make_launch_context();
    (*ker)(launch_ctx);
    std::cout << prog_.fetch_result<int>(0) << " : ans\n";
    exit(0);
}

void test_exptr() {
    auto prog_ = Program(arch_from_name("x64"));
    //prog_.thread_pool = std::make_unique<ThreadPool>(1);
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    prog_.config = config_print_ir;  // print ir
    prog_.materialize_layout();

    IRBuilder builder;
    auto *ten = builder.get_float32(10);
    auto *index = builder.get_int32(1);
    auto *arg = builder.create_arg_load(0, PrimitiveType::f32, true);
    auto *ptr = builder.insert(std::make_unique<ExternalPtrStmt>(arg, std::vector<Stmt*>(1, index)));
    builder.insert(std::make_unique<GlobalStoreStmt>(ptr, ten));

    auto block = builder.extract_ir();
    auto ker = std::make_unique<Kernel>(prog_, std::move(block));
    ker->insert_arg(PrimitiveType::gen, true);
    auto launch_ctx = ker->make_launch_context();
    float tmp[10] = {0, 1, 2, 3};
    launch_ctx.set_arg_external_array(0, uint64(tmp), 4);
    (*ker)(launch_ctx);
    //std::cout << prog_.fetch_result<int>(0) << " : ans\n";
    std::cout << tmp[0] << " : tmp[0]\n";
    std::cout << tmp[1] << " : tmp[1]\n";
    std::cout << tmp[2] << " : tmp[2]\n";
    std::cout << tmp[3] << " : tmp[3]\n";
    //std::cout << int(tmp[0]) << " : int tmp[0]\n";
    exit(0);
}

class GradInfoChiPrimal final : public SNode::GradInfoProvider {
public:
    SNode &grad;
    explicit GradInfoChiPrimal(SNode &_grad): grad(_grad) {}
    bool is_primal() const override { return true; }
    SNode *grad_snode() const override { return &grad; }
};

class GradInfoChiAdjoint final : public SNode::GradInfoProvider {
public:
    explicit GradInfoChiAdjoint() {}
    bool is_primal() const override { return false; }
    SNode *grad_snode() const override { return nullptr; }
};

void test_auto_diff() {
    auto prog_ = Program(arch_from_name("x64"));
    //prog_.thread_pool = std::make_unique<ThreadPool>(1);
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    prog_.config = config_print_ir;  // print ir
    auto root = prog_.snode_root.get();

    auto grad_snode = [&](Index index, int size) -> SNode& {
        auto &snode = root->dense(index, size).insert_children(SNodeType::place);
        snode.dt = PrimitiveType::f32;
        snode.grad_info = std::make_unique<GradInfoChiPrimal>(root->dense(index, size).insert_children(SNodeType::place));
        snode.get_grad()->dt = PrimitiveType::f32;
        snode.get_grad()->grad_info = std::make_unique<GradInfoChiAdjoint>();
        return snode;
    };
    auto &energy = grad_snode(0, 1);
    auto &x = grad_snode(0, 2);
    prog_.materialize_layout();

    {  // init
      IRBuilder builder;
      builder.create_global_store(
          builder.create_global_ptr(
              &x, std::vector<Stmt *>(1, builder.get_int32(0))),
          builder.get_int32(7));
      builder.create_global_store(
          builder.create_global_ptr(
              &x, std::vector<Stmt *>(1, builder.get_int32(1))),
          builder.get_int32(23));
      builder.create_global_store(
          builder.create_global_ptr(
              energy.get_grad(), std::vector<Stmt *>(1, builder.get_int32(1))),
          builder.get_int32(1));

      auto ker = std::make_unique<Kernel>(prog_, builder.extract_ir());
      auto launch_ctx = ker->make_launch_context();
      (*ker)(launch_ctx);
    }
    {  // grad
      IRBuilder builder;
      auto *dst = builder.create_global_ptr(
          &energy, std::vector<Stmt *>(1, builder.get_int32(0)));
      auto *val_0 = builder.create_global_load(builder.create_global_ptr(
          &x, std::vector<Stmt *>(1, builder.get_int32(0))));
      auto *val_1 = builder.create_global_load(builder.create_global_ptr(
          &x, std::vector<Stmt *>(1, builder.get_int32(1))));
      auto *val = builder.create_mul(val_0, val_1);
      builder.insert(
          std::make_unique<AtomicOpStmt>(AtomicOpType::add, dst, val));

      auto block = builder.extract_ir();
      auto ker =
          std::make_unique<Kernel>(prog_, std::move(block), "auto_diff", true);
      // auto ker = std::make_unique<Kernel>(prog_, std::move(block));
      auto launch_ctx = ker->make_launch_context();
      (*ker)(launch_ctx);
    }
    {  // output
      IRBuilder builder;
      builder.create_return(
          builder.create_global_load(builder.create_global_ptr(
              x.get_grad(), std::vector<Stmt *>(1, builder.get_int32(0)))));

      auto ker = std::make_unique<Kernel>(prog_, builder.extract_ir());
      auto launch_ctx = ker->make_launch_context();
      (*ker)(launch_ctx);
      printf("ans = %lf\n", prog_.fetch_result<float>(0));
    }
    exit(0);
}
}  // namespace lang
}  // namespace taichi

namespace ra {
std::unique_ptr<taichi::lang::Kernel> ker_init, ker, ker_grad, ker_output;
std::unique_ptr<taichi::GUI> gui_ptr;
std::set<std::string> names, reserved_funcs;
const int W = 512, H = 512, thread_num = 8;
taichi::real *args = nullptr;
std::unique_ptr<taichi::lang::Program> program_ptr;

void init() {
    ker_init.reset(nullptr);
    ker.reset(nullptr);
    ker_grad.reset(nullptr);
    ker_output.reset(nullptr);
    gui_ptr.reset(nullptr);
    program_ptr.reset(nullptr);
    names.clear();
    reserved_funcs.clear();
    delete [] args;
}

void screenshot() {
    gui_ptr->screenshot();
}

void run_string(std::string plain) {
    //taichi::lang::test_snode();
    //taichi::lang::test_exptr();
    //taichi::lang::test_auto_diff();
    /*std::fstream fs(argv[1], std::fstream::in);
    std::string plain;
    for (char ch; fs.get(ch); plain += ch);
    fs.close();*/
    reserved_funcs.insert("cos");
    reserved_funcs.insert("sin");
    auto &&lex = get_lexical(plain, names, reserved_funcs);
    names.erase("x");
    names.erase("y");
    //printf("plain: %d, lex: %d\n", plain.size(), lex.size());
    auto ast_root = std::unique_ptr<Expr>(get_ast(lex));
    print(ast_root);
    //fs.open(argv[2], std::fstream::out);
    //init_IR(fs, names);
    //int cnt = 0;
    //dummy_IR(root, fs, cnt);
    //fs << "\treturn ssa[" + std::to_string(cnt - 1) + "]";
    //fs << "\treturn ssa" + std::to_string(cnt - 1) + "";
    //delete program_ptr.get();
    program_ptr = std::make_unique<taichi::lang::Program>(taichi::lang::arch_from_name("cuda"));
    std::cerr << "OK\n";
    {
        using namespace taichi::lang;
        auto &program = *program_ptr;
        CompileConfig config_print_ir;
        config_print_ir.print_ir = true;
        program.config = config_print_ir;  // print ir
        auto root = program.snode_root.get();
        std::vector<Index> index_dense = {0, 1};
        std::vector<int> size_dense = {W, H};
        auto grad_snode = [&](std::vector<Index> index, std::vector<int> size) -> SNode& {
            auto &snode = root->dense(index, size).insert_children(SNodeType::place);
            snode.dt = PrimitiveType::f32;
            snode.grad_info = std::make_unique<GradInfoChiPrimal>(
                root->dense(index, size).insert_children(SNodeType::place));
            snode.get_grad()->dt = PrimitiveType::f32;
            snode.get_grad()->grad_info = std::make_unique<GradInfoChiAdjoint>();
            return snode;
        };
        auto &place = grad_snode(index_dense, size_dense);
        auto &matrix_x = grad_snode(index_dense, size_dense);
        auto &matrix_y = grad_snode(index_dense, size_dense);
        program.materialize_layout();
        {  // init
          IRBuilder builder;
          auto *zero = builder.get_int32(0);
          auto *width = builder.get_int32(W);
          auto *height = builder.get_int32(H);
          auto *loopx = builder.create_range_for(zero, width, 1, 0, thread_num);
          {
            auto guard_x = builder.get_loop_guard(loopx);
            auto *index_x = builder.get_loop_index(loopx, 0);
            auto *x = builder.create_div(
                builder.create_cast(index_x, PrimitiveType::f32),
                builder.create_cast(width, PrimitiveType::f32));
            auto *loopy =
                builder.create_range_for(zero, height, 1, 0, thread_num);
            { 
                auto guard_y = builder.get_loop_guard(loopy); 
                auto *index_y = builder.get_loop_index(loopy, 0);
                auto *y = builder.create_div(builder.create_cast(index_y, PrimitiveType::f32),
                    builder.create_cast(height, PrimitiveType::f32));
                std::vector<Stmt*> indices = {index_x, index_y};
                builder.create_global_store(
                    builder.create_global_ptr(place.get_grad(), indices),
                    builder.get_float32(1));
                builder.create_global_store(
                    builder.create_global_ptr(matrix_x.get_grad(), indices),
                    builder.get_float32(0));
                builder.create_global_store(
                    builder.create_global_ptr(matrix_y.get_grad(), indices),
                    builder.get_float32(0));
                builder.create_global_store(
                    builder.create_global_ptr(&matrix_x, indices), x);
                builder.create_global_store(
                    builder.create_global_ptr(&matrix_y, indices), y);
            }
          }
          ker_init = std::make_unique<Kernel>(program, builder.extract_ir());
        }
        std::cerr << "after init\n";

        auto get_main_ir = [&](bool grad) -> auto {
            IRBuilder builder;
            std::map<std::string, Stmt*> para;
            {
                int i = 0;
                for (auto name : names) {
                    para[name] = builder.create_arg_load(i + 1, PrimitiveType::f32, false); // arg 0 is canvas
                    i++;
                }
            }
            auto *zero = builder.get_int32(0);
            auto *width = builder.get_int32(W);
            auto *height = builder.get_int32(H);
            auto *loopx = builder.create_range_for(zero, width, 1, 0, thread_num);
            {
                auto guard_x = builder.get_loop_guard(loopx);
                auto *index_x = builder.get_loop_index(loopx, 0);
                auto *loopy = builder.create_range_for(zero, height, 1, 0, thread_num);
                {
                    auto guard_y = builder.get_loop_guard(loopy);
                    auto *index_y = builder.get_loop_index(loopy, 0);
                    std::vector<Stmt*> indices = {index_x, index_y};
                    para["x"] = builder.create_global_load(builder.create_global_ptr(&matrix_x, indices));
                    para["y"] = builder.create_global_load(builder.create_global_ptr(&matrix_y, indices));
                    auto *ans = build_IR(ast_root, builder, para);
                    { // snode
                        builder.create_global_store(builder.create_global_ptr(&place, indices), ans);
                    }
                    /*{ // external
                        auto *index = builder.create_add(builder.create_mul(index_x, height), index_y);
                        builder.create_global_store(
                            builder.create_external_ptr(
                                builder.create_arg_load(0, PrimitiveType::f32,
                                                        true),
                                std::vector<Stmt *>(1, index)),
                            ans);
                    }*/
                }
            }
            auto ker = new Kernel(program, builder.extract_ir(), "main", grad);
            ker->insert_arg(PrimitiveType::gen, true);
            for (auto name : names) {
                ker->insert_arg(PrimitiveType::f32, false);
            }
            return ker;
        };
        std::cerr << "after main\n";
        ker = std::unique_ptr<Kernel>(get_main_ir(false));
        ker_grad = std::unique_ptr<Kernel>(get_main_ir(true));
        {  // output
          IRBuilder builder;
          auto *zero = builder.get_int32(0);
          auto *width = builder.get_int32(W);
          auto *height = builder.get_int32(H);
          auto *loopx = builder.create_range_for(zero, width, 1, 0, thread_num);
          {
            auto guard_x = builder.get_loop_guard(loopx);
            auto *index_x = builder.get_loop_index(loopx, 0);
            auto *loopy =
                builder.create_range_for(zero, height, 1, 0, thread_num);
            { 
                auto guard_y = builder.get_loop_guard(loopy); 
                auto *index_y = builder.get_loop_index(loopy, 0);
                std::vector<Stmt*> indices = {index_x, index_y};
                auto *index = builder.create_add(builder.create_mul(index_x, height), index_y);
                //auto ans = builder.create_global_load(builder.create_global_ptr(matrix_x.get_grad(), indices));
                auto ans = builder.create_global_load(builder.create_global_ptr(&place, indices));
                builder.create_global_store(
                    builder.create_external_ptr(
                        builder.create_arg_load(0, PrimitiveType::f32, true),
                        std::vector<Stmt *>(1, index)),
                    ans);
                auto grad_x = builder.create_global_load(builder.create_global_ptr(matrix_x.get_grad(), indices));
                builder.create_global_store(
                    builder.create_external_ptr(
                        builder.create_arg_load(1, PrimitiveType::f32, true),
                        std::vector<Stmt *>(1, index)),
                    grad_x);
                auto grad_y = builder.create_global_load(builder.create_global_ptr(matrix_y.get_grad(), indices));
                builder.create_global_store(
                    builder.create_external_ptr(
                        builder.create_arg_load(2, PrimitiveType::f32, true),
                        std::vector<Stmt *>(1, index)),
                    grad_y);
            }
          }
          ker_output = std::make_unique<Kernel>(program, builder.extract_ir());
          for (int i = 0; i < 3; i++) ker_output->insert_arg(PrimitiveType::gen, true);
        }
        //delete gui_ptr.get();
        auto ptr = std::make_unique<taichi::GUI>("GUI Test", 720, H, true, false, 0, false, false);
        ptr->update();
        gui_ptr = std::make_unique<taichi::GUI>("GUI Test", 720, H, true, false, 0, false, false);
        auto &gui = *gui_ptr;
        args = new taichi::real[names.size()];
        memset(args, 0, sizeof(taichi::real) * names.size());
        {
            int i = 0;
            for (auto name : names) {
                gui.slider(name, args[i], taichi::real(0), taichi::real(1));
                i++;
            }
        }
        //gui.button("sc", screenshot);
    }
}

void run_frame() {
    auto &gui = *gui_ptr;
    auto &canvas = *gui.canvas;
    static float ans[W][H], grad_x[W][H], grad_y[W][H];
    auto launch_ctx = ker->make_launch_context();
    auto ctx_init = ker_init->make_launch_context();
    auto ctx_grad = ker_grad->make_launch_context();
    auto ctx_output = ker_output->make_launch_context();
    (*ker_init)(ctx_init);
    launch_ctx.set_arg_external_array(0, taichi::uint64(ans), 0);
    {
        int i = 0;
        for (auto name : names) {
            launch_ctx.set_arg_float(i + 1, args[i]);
            i++;
        }
    } 
    (*ker)(launch_ctx);
    ctx_grad.set_arg_external_array(0, taichi::uint64(ans), 0);
    {
        int i = 0;
        for (auto name : names) {
            ctx_grad.set_arg_float(i + 1, args[i]);
            i++;
        }
    }
    (*ker_grad)(ctx_grad);
    ctx_output.set_arg_external_array(0, taichi::uint64(ans), 0);
    ctx_output.set_arg_external_array(1, taichi::uint64(grad_x), 0);
    ctx_output.set_arg_external_array(2, taichi::uint64(grad_y), 0);
    (*ker_output)(ctx_output);
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            float r = 0, g = 0, b = 0;
            double tmp = ans[i][j] * 256;
            //double tmp = (grad_y[i][j] - grad_x[i][j]) * 256;
            if (tmp < 32) b = 128 + tmp * 4;
            else if (tmp < 96) b = 255;
            else if (tmp < 159) b = 254 - (tmp - 96) * 4;
            if (tmp < 33) g = 0;
            else if (tmp < 95) g = 4 * (tmp - 32);
            else if (tmp < 160) g = 255;
            else if (tmp < 223) g = 252 - 4 * (tmp - 160);
            if (tmp > 224) r = 252 - 4 * (tmp - 224);
            else if (tmp > 158) r = 255;
            else if (tmp > 96) r = 4 * (tmp - 97);
            r = r / 255;
            g = g / 255;
            b = b / 255;
            //std::array<taichi::real, 4> tmp{0, ans[i][j] * (1 - blue), ans[i][j] * blue, 255};
            std::array<taichi::real, 4> color{r, g, b, 1};
            //std::array<taichi::real, 4> tmp{0, (float)(i + j) / W, 0, 255};
            canvas.img[i][j] = taichi::Vector4(color);
        }
    }
    for (int i = 0; i < W; i += 32) {
        for (int j = 0; j < H; j += 32) {
            using namespace taichi;
            float k = 3, a = 0.2, w = 2, arrow_k = 0.8;
            auto off = Vector2(grad_x[i][j] * k, grad_y[i][j] * k);
            auto left = Vector2(off[0] * std::cos(a) - off[1] * std::sin(a),
                        off[0] * std::sin(a) + off[1] * std::cos(a)) * arrow_k;
            auto right = Vector2(off[0] * std::cos(-a) - off[1] * std::sin(-a),
                        off[0] * std::sin(-a) + off[1] * std::cos(-a)) * arrow_k;
            auto src = Vector2(i, j);
            canvas.path(src, src + off)
                .close().color(0, 0, 0).width(w).finish();
            canvas.path(src + off, src + left)
                .close().color(0, 0, 0).width(w).finish();
            canvas.path(src + off, src + right)
                .close().color(0, 0, 0).width(w).finish();
        }
    }
    gui.update();
    for (;gui.has_key_event();) {
        auto eve = gui.get_key_event_head();
        gui.pop_key_event_head();
        if (eve.type == taichi::GUI::KeyEvent::Type::press && eve.key == "s") {
            screenshot();
        }
    }
}
}

int main() {
    using namespace std;
    string input = "";
    for (char ch; ~(ch = getchar());) {
        input += ch;
    }
    ra::init();
    ra::run_string(input);
    for (;;) ra::run_frame();
    return 0;
}
