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

void game_of_life() {
    using namespace std;
    using namespace taichi;
    using namespace lang;
    auto program = Program(arch_from_name("x64"));
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    program.config = config_print_ir;  // print ir
    auto root = program.snode_root.get();
    auto &dense = root->dense(Index(0), 10);
    auto &alive = dense.insert_children(SNodeType::place);
    alive.dt = PrimitiveType::i32;
    auto &next = dense.insert_children(SNodeType::place);
    next.dt = PrimitiveType::i32;
    program.materialize_layout();

    const int N = 512;
    std::unique_ptr<Kernel> step_kernel;
    {
        IRBuilder builder;
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
                Stmt *sum = nullptr, *sum_ex = nullptr;
                for (int i = -1; i < 2; i++) {
                    for (int j = -1; j < 2; j++) {
                        vector<Stmt *> indices = {
                            builder.create_add(builder.get_int32(i), x),
                            builder.create_add(builder.get_int32(j), y)};
                        auto *load = builder.create_global_load(builder.create_global_ptr(&alive, indices));
                        if (sum) sum = builder.create_add(sum, load);
                        else sum = load;
                        if (!i && !j) sum_ex = builder.create_sub(sum, load);
                    }
                }
                auto ans = builder.create_or(
                    builder.create_cmp_eq(sum, builder.get_int32(3)),
                    builder.create_cmp_eq(sum_ex, builder.get_int32(3)));
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(builder.create_global_ptr(&next, indices), ans);
            }
        }
        step_kernel = make_unique<Kernel>(program, builder.extract_ir());
    }

    std::unique_ptr<Kernel> swap_kernel;
    {
        IRBuilder builder;
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
        swap_kernel = make_unique<Kernel>(program, builder.extract_ir());
    }
}

int main() {
    game_of_life();
    return 0;
}