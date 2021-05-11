#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <stack>
#include <set>
#include <cassert>

#include "taichi/ir/ir_builder.h"
#include "taichi/ir/statements.h"
#include "taichi/program/program.h"
#include "taichi/ir/ir.h"

#include "taichi/gui/gui.h"

class BuilderHelper {
public:
    static taichi::lang::IRBuilder *builder;
    taichi::lang::Stmt *stmt;

    BuilderHelper(taichi::lang::Stmt *_stmt): stmt(_stmt) {}

    BuilderHelper operator + (taichi::lang::Stmt *other) {
        return builder->create_add(stmt, other);
    }

    BuilderHelper operator - (taichi::lang::Stmt *other) {
        return builder->create_sub(stmt, other);
    }

    BuilderHelper operator * (taichi::lang::Stmt *other) {
        return builder->create_mul(stmt, other);
    }

    BuilderHelper operator / (taichi::lang::Stmt *other) {
        return builder->create_div(stmt, other);
    }

    BuilderHelper operator % (taichi::lang::Stmt *other) {
        return builder->create_mod(stmt, other);
    }

    BuilderHelper operator > (taichi::lang::Stmt *other) {
        return builder->create_cmp_gt(stmt, other);
    }

    BuilderHelper operator | (taichi::lang::Stmt *other) {
        return builder->create_or(stmt, other);
    }

    BuilderHelper operator == (taichi::lang::Stmt *other) {
        return builder->create_cmp_eq(stmt, other);
    }

    BuilderHelper operator & (taichi::lang::Stmt *other) {
        return builder->create_and(stmt, other);
    }

    operator taichi::lang::Stmt*() {
        return stmt;
    }
};
taichi::lang::IRBuilder *BuilderHelper::builder(nullptr);

class BuilderHelperGuard {
public:
    BuilderHelperGuard(taichi::lang::IRBuilder &builder) {
        BuilderHelper::builder = &builder;
    }

    ~BuilderHelperGuard() {
        BuilderHelper::builder = nullptr;
    }
};

void game_of_life() {
    using namespace std;
    using namespace taichi;
    using namespace lang;
    auto program = Program(arch_from_name("x64"));
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    program.config = config_print_ir;  // print ir
    auto root = program.snode_root.get();

    const int N = 512;
    std::vector<Index> index_dense = {0, 1};
    std::vector<int> size_dense = {N, N};
    auto &alive = root->dense(index_dense, size_dense).insert_children(SNodeType::place);
    alive.dt = PrimitiveType::i32;
    auto &next = root->dense(index_dense, size_dense).insert_children(SNodeType::place);
    next.dt = PrimitiveType::i32;
    program.materialize_layout();

    typedef BuilderHelper BH;

    std::unique_ptr<Kernel> kernel_init;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *left = builder.get_int32(1);
        auto *right = builder.get_int32(N - 1);
        auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop_x);
            auto *x = builder.get_loop_index(loop_x, 0);
            auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_y);
                auto *y = builder.get_loop_index(loop_y, 0);
                Stmt *cnt = BH(x) * builder.get_int32(N) + y;
                auto *mod = builder.get_int32(10007);
                Stmt *ans = (BH(cnt) * builder.get_int32(9287) % mod + mod) % builder.get_int32(2);
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(builder.create_global_ptr(&alive, indices), ans);
            }
        }
        kernel_init = make_unique<Kernel>(program, builder.extract_ir());
    }
    
    std::unique_ptr<Kernel> kernel_step;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *left = builder.get_int32(1);
        auto *right = builder.get_int32(N - 1);
        auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop_x);
            auto *x = builder.get_loop_index(loop_x, 0);
            auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_y);
                auto *y = builder.get_loop_index(loop_y, 0);
                Stmt *sum = nullptr, *self = nullptr;
                for (int i = -1; i < 2; i++) {
                    for (int j = -1; j < 2; j++) {
                        vector<Stmt *> indices = {
                            builder.create_add(builder.get_int32(i), x),
                            builder.create_add(builder.get_int32(j), y)};
                        auto *load = builder.create_global_load(builder.create_global_ptr(&alive, indices));
                        if (sum) sum = builder.create_add(sum, load);
                        else sum = load;
                        if (!i && !j) self = load;
                    }
                }
                auto ans = (BH(sum) == builder.get_int32(3) | (BH(sum) - self) == builder.get_int32(3)) & builder.get_int32(1);
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(builder.create_global_ptr(&next, indices), ans);
            }
        }
        kernel_step = make_unique<Kernel>(program, builder.extract_ir());
    }

    std::unique_ptr<Kernel> kernel_swap;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *left = builder.get_int32(1);
        auto *right = builder.get_int32(N - 1);
        auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop_x);
            auto *x = builder.get_loop_index(loop_x, 0);
            auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_y);
                auto *y = builder.get_loop_index(loop_y, 0);
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(
                    builder.create_global_ptr(&alive, indices),
                    builder.create_global_load(
                        builder.create_global_ptr(&next, indices)));
            }
        }
        kernel_swap = make_unique<Kernel>(program, builder.extract_ir());
    }

    std::unique_ptr<Kernel> kernel_gui;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *left = builder.get_int32(1);
        auto *right = builder.get_int32(N - 1);
        auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop_x);
            auto *x = builder.get_loop_index(loop_x, 0);
            auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_y);
                auto *y = builder.get_loop_index(loop_y, 0);
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(
                    builder.create_external_ptr(
                        builder.create_arg_load(0, PrimitiveType::i32, true),
                        vector<Stmt *>(1, BH(x) * builder.get_int32(N) + y)),
                    builder.create_global_load(
                        builder.create_global_ptr(&alive, indices)));
            }
        }
        kernel_gui = make_unique<Kernel>(program, builder.extract_ir());
        kernel_gui->insert_arg(PrimitiveType::gen, true);
    }
    auto gui = GUI("cell language", N, N, true, false, 0, false, false);
    auto &canvas = *gui.canvas;
    static int ans[N][N];
    auto ctx_init = kernel_init->make_launch_context();
    auto ctx_step = kernel_step->make_launch_context();
    auto ctx_swap = kernel_swap->make_launch_context();
    auto ctx_gui = kernel_gui->make_launch_context();
    ctx_gui.set_arg_external_array(0, taichi::uint64(ans), 0);
    (*kernel_init)(ctx_init);
    for (int frame = 0;; frame++) {
        gui.update();
        if (frame % 60) continue;
        (*kernel_step)(ctx_step);
        if (frame > 60) (*kernel_swap)(ctx_swap);
        (*kernel_gui)(ctx_gui);
        long long sum = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                float k = ans[i][j];
                sum += ans[i][j];
                std::array<taichi::real, 4> color{k, 0, 0, 1};
                canvas.img[i][j] = taichi::Vector4(color);
            }
        }
        cerr << sum << ": sum\n";
    }
}

int main() {
    game_of_life();
    return 0;
}