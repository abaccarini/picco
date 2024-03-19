/*
   PICCO: A General Purpose Compiler for Private Distributed Computation
   ** Copyright (C) from 2013 PICCO Team
   ** Department of Computer Science and Engineering, University of Notre Dame

   PICCO is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PICCO is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with PICCO. If not, see <http://www.gnu.org/licenses/>.
*/

/* This file is modified from OMPi */

/* ast_show.c -- prints out the tree; makes it look good, too. */

#include "ast_show.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static symbol _curfile = NULL;
static int indlev = 0; /* Indentation level */
int ind = 0;
int priv_if_index = 0;
int delete_tmp_array[3 * 1000];
int buffer_size = 100;
condstack if_top = NULL;
branchnode if_tree = NULL;
batch_condition_stack bcs = NULL;
batch_statement_stack bss = NULL;
batch_private_index_stack bpis = NULL;
control_sequence_stack batch_stack = NULL;
control_sequence_stack private_selection_stack = NULL;
struct_node_stack sns = NULL;
FILE *output = NULL;
FILE *global_stream = NULL; /* The stream that will hold global private variables*/
int gf = 0; /* A global flag that will be updated based on the tree->gflag. 1-Global, 0-Regular*/
int is_priv = 0; /* A global flag that will be updated based on the varible. 1-Private, 0-Public*/
is_smc_io = 0; /* Used to print the name and values of variables to two streams 1-print for smc_input/O, 0-Not*/
is_smc_set = 0; /* Used to print the name and values of variables to two streams 1-print for smc_set, 0-Not*/
int global_batch_flag = 0;
int declared = 0;
int enterfunc = 0;
int is_return_void = 0;
str arg_str;

static void indent() {
    int i;
    if (indlev > 0)
        for (i = 2 * indlev; i > 0; i--)
            putc(' ', output);
}

static void indent_global_stream() { 
    int i;
    if (indlev > 0)
        for (i = 2 * indlev; i > 0; i--)
            putc(' ', global_stream);
}

void ast_clear_all_tmp_variables(branchnode current) {
    ast_free_memory_for_local_variables(current->tablelist);
    if (tmp_index >= 1)
        ast_tmp_clear_show("_picco_tmp", 1, tmp_index);
    if (tmp_float_index >= 1)
        ast_float_tmp_clear_show("_picco_ftmp", 1, tmp_float_index);

    if (is_priv_int_index_appear || is_priv_float_index_appear)
        ast_tmp_clear_show("_picco_priv_ind", 1, 3);
    if (is_priv_int_index_appear || is_priv_int_ptr_appear)
        ast_tmp_clear_show("_picco_priv_tmp", 1, 2);
    if (is_priv_float_index_appear || is_priv_float_ptr_appear) {
        indent();
        fprintf(output, "for(int i = 0; i < 4; i++){\n");
        indlev++;
        indent();
        fprintf(output, "ss_clear(_picco_priv_ftmp1[i]);\n");
        indent();
        fprintf(output, "ss_clear(_picco_priv_ftmp2[i]);\n");
        indlev--;
        fprintf(output, "}\n");
        indent();
        fprintf(output, "free(_picco_priv_ftmp1);\n");
        indent();
        fprintf(output, "free(_picco_priv_ftmp2);\n");
    }
    if (is_priv_int_ptr_appear) {
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_int_ptr1);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_int_ptr2);\n");
    }
    if (is_priv_float_ptr_appear) {
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_float_ptr1);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_float_ptr2);\n");
    }

    if (is_priv_int_struct_field_appear) {
        indent();
        fprintf(output, "ss_clear(_picco_str_field_tmp_int1);\n");
        indent();
        fprintf(output, "ss_clear(_picco_str_field_tmp_int2);\n");
        indent();
        fprintf(output, "ss_clear(_picco_str_field_tmp_int3);\n");
    }

    if (is_priv_int_ptr_struct_field_appear) {
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_int_ptr1);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_int_ptr2);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_int_ptr3);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_int_ptr);\n");
    }

    if (is_priv_float_struct_field_appear) {
        indent();
        fprintf(output, "for(int _picco_i = 0; _picco_i < 4; _picco_i++){\n");
        indlev++;
        indent();
        fprintf(output, "ss_clear(_picco_str_field_tmp_float1[_picco_i]);\n");
        indent();
        fprintf(output, "ss_clear(_picco_str_field_tmp_float2[_picco_i]);\n");
        indent();
        fprintf(output, "ss_clear(_picco_str_field_tmp_float3[_picco_i]);\n");
        indlev--;
        fprintf(output, "}\n");
        indent();
        fprintf(output, "free(_picco_str_field_tmp_float1);\n");
        indent();
        fprintf(output, "free(_picco_str_field_tmp_float2);\n");
        indent();
        fprintf(output, "free(_picco_str_field_tmp_float3);\n");
    }

    if (is_priv_float_ptr_struct_field_appear) {
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_float_ptr1);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_float_ptr2);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_float_ptr3);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_float_ptr);\n");
    }

    if (is_priv_struct_ptr_struct_field_appear) {
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_struct_ptr1);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_struct_ptr2);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_struct_ptr3);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_struct_ptr4);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_struct_ptr5);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_str_field_tmp_struct_ptr6);\n");
        indent();
        fprintf(output, "__s->smc_free_ptr(&_picco_tmp_struct_ptr);\n");
    }
    //  ast_tmp_clear_show("cond", 1, 3);
}

void ast_stmt_jump_show(aststmt tree, branchnode current) {
    switch (tree->subtype) {
    case SBREAK:
        fprintf(output, "break;");
        break;
    case SCONTINUE:
        fprintf(output, "continue;");
        break;
    case SRETURN:
        ast_clear_all_tmp_variables(current);
        // if return value is public
        // if(tree->u.expr == NULL || tree->u.expr->flag == PUB)
        if (tree->u.expr == NULL || tree->u.expr->flag == PUB || (tree->u.expr->type == CASTEXPR)) {
            indent();
            fprintf(output, "return");
            if (tree->u.expr != NULL) {
                fprintf(output, " (");
                ast_expr_show(tree->u.expr);
                fprintf(output, ")");
            }
        }
        // return type is private
        // assume for now, the return value is always an identifier
        else {
            indent();
            fprintf(output, "ss_set(rvar, %s)", tree->u.expr->u.sym->name);
        }
        fprintf(output, ";");
        // if return value is private
        break;
    case SGOTO:
        fprintf(output, "goto %s;", tree->u.label->name);
        break;
    default:
        fprintf(stderr, "[ast_stmt_jump_show]: b u g !!\n");
    }
    fprintf(output, "\n");
}

void ast_stmt_batch_show(aststmt tree, branchnode current) {
    int batch_index = 0;
    int statement_index = 0;
    int narray_element_index = 0;
    int private_index = 0;
    int private_selection_index = 0;
    batch_stack = control_sequence_stack_new();
    private_selection_stack = control_sequence_stack_new();
    ast_batch_iter_tree(tree, &batch_index, &statement_index, &private_selection_index);
    fprintf(output, "{\n");
    indlev++;
    indent();
    ast_batch_declare_counter(tree, "_picco_batch_counter", batch_index, current);
    ast_batch_declare_counter(tree, "_picco_ind", batch_index, current);
    batch_index = 0;
    ast_batch_compute_counter(tree, &batch_index, current);
    ast_batch_allocate_counter();
    // create an array for holding the variables in the expression
    batch_index = 0;
    ast_batch_declare_array_for_narrayelement(tree, &narray_element_index, &private_index, &batch_index);
    batch_index = 0;
    statement_index = 0;
    narray_element_index = 0;
    private_index = 0;

    ast_batch_compute_index(tree, &batch_index, &statement_index, &narray_element_index, delete_tmp_array, &private_index);

    // re-initialize all control variables.
    batch_index = 0;
    statement_index = 0;
    private_selection_index = 0;
    narray_element_index = 0;
    private_index = 0;
    ast_batch_compute_stmt(tree, &batch_index, &statement_index, &private_selection_index, &narray_element_index, &private_index, current);
    ast_batch_clear_counter(private_selection_index, narray_element_index, delete_tmp_array);
    indlev--;
    indent();
    fprintf(output, "}\n");
    batch_statement_popAll(bss);
}

void ast_batch_compute_stmt(aststmt tree, int *batch_index, int *statement_index, int *private_selection_index, int *narray_element_index, int *private_index, branchnode current) {
    int tmp_index = 0;
    switch (tree->type) {
    case COMPOUND:
        ast_batch_compute_stmt(tree->body, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
        break;
    case STATEMENTLIST:
        ast_batch_compute_stmt(tree->u.next, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
        ast_batch_compute_stmt(tree->body, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
        break;
    case BATCH:
        (*batch_index)++;
        control_sequence_push(*batch_index, batch_stack);
        ast_batch_compute_stmt(tree->body, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
        control_sequence_pop(batch_stack);
        break;
    case ITERATION:
        ast_batch_compute_stmt(tree->body, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
        break;
    case SELECTION:
        if (tree->u.selection.cond->flag == PUB) {
            (*batch_index)++;
            control_sequence_push(*batch_index, batch_stack);
            ast_batch_compute_stmt(tree->body, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
            control_sequence_pop(batch_stack);
            if (tree->u.selection.elsebody) {
                (*batch_index)++;
                control_sequence_push(*batch_index, batch_stack);
                ast_batch_compute_stmt(tree->u.selection.elsebody, batch_index, statement_index, private_selection_index, narray_element_index, private_index, current);
                control_sequence_pop(batch_stack);
            }
        }

        if (tree->u.selection.cond->flag == PRI) {
            int batch_value = *batch_index;
            (*private_selection_index)++;
            tmp_index = *private_selection_index;
            (*statement_index)++;
            control_sequence_push(*private_selection_index, private_selection_stack);
            private_selection_stack->head->batch_index = batch_value;
            branchnode left = NULL;
            left = if_branchnode_insert(current, NULL, 0, *private_selection_index, -1, 0, 0);
            ast_batch_print_stmt(Expression(tree->u.selection.cond), batch_stack->head->index, *statement_index, narray_element_index, private_index, left);
            ast_batch_compute_stmt(tree->body, batch_index, statement_index, private_selection_index, narray_element_index, private_index, left);
            if_branchtree_remove(current, left);
            control_sequence_pop(private_selection_stack);

            if (tree->u.selection.elsebody) {
                (*private_selection_index)++;
                (*statement_index)++;
                control_sequence_push(*private_selection_index, private_selection_stack);
                private_selection_stack->head->batch_index = batch_value;
                branchnode right = NULL;
                right = if_branchnode_insert(current, NULL, 1, *private_selection_index, tmp_index, 0, 0);
                ast_batch_print_stmt(Expression(tree->u.selection.cond), batch_stack->head->index, *statement_index, narray_element_index, private_index, right);
                ast_batch_compute_stmt(tree->u.selection.elsebody, batch_index, statement_index, private_selection_index, narray_element_index, private_index, right);
                if_branchtree_remove(current, right);
                control_sequence_pop(private_selection_stack);
            }
        }
        break;
    case EXPRESSION:
        switch (tree->u.expr->type) {
        case ASS:
            (*statement_index)++;
            ast_batch_print_stmt(tree, batch_stack->head->index, *statement_index, narray_element_index, private_index, current);
            break;
        }
        break;
    default:
        break;
    }
}

void ast_print_private_indexed_stmt_write(astexpr tree, char *input_index, char *output_result, int batch_index, int opid, char *buf1, char *buf2, char *buf3) {
    char *type = (char *)malloc(sizeof(char) * buffer_size);
    str arrayname = Str("");
    str arraysize1 = Str("");
    str arraysize2 = Str("");

    if (opid != 0) {
        printf("Please use simple assignment operators '='.\n");
        exit(0);
    }

    if (tree->ftype == 0)
        sprintf(type, "int");
    else if (tree->ftype == 1)
        sprintf(type, "float");

    // one-dimension
    if (tree->left->type == IDENT) {
        ast_expr_print(arrayname, tree->left);
        ast_expr_print(arraysize1, tree->left->arraysize);
        indent();
        if (batch_index != -1)
            fprintf(output, "__s->smc_privindex_write(%s, %s, %d, %d, %s, %s, _picco_batch_counter%d, %s, %s, %s, \"%s\", %d)", input_index, str_string(arrayname), tree->size, tree->sizeexp, output_result, str_string(arraysize1), batch_index, buf1, buf2, buf3, type, tree->thread_id);
        else
            fprintf(output, "__s->smc_privindex_write(%s, %s, %d, %d, %s, %s, %s, NULL, -1, \"%s\", %d)", input_index, str_string(arrayname), tree->size, tree->sizeexp, output_result, str_string(arraysize1), buf1, type, tree->thread_id);
    }
    // two-dimension
    else {
        if (tree->left->type == ARRAYIDX) {
            ast_expr_print(arrayname, tree->left->left);
            ast_expr_print(arraysize1, tree->left->arraysize);
            ast_expr_print(arraysize2, tree->left->left->arraysize);
            indent();
            if (batch_index != -1)
                fprintf(output, "__s->smc_privindex_write(%s, %s, %d, %d, %s, %s, %s, _picco_batch_counter%d, %s, %s, %s, \"%s\", %d)", input_index, str_string(arrayname), tree->size, tree->sizeexp, output_result, str_string(arraysize1), str_string(arraysize2), batch_index, buf1, buf2, buf3, type, tree->thread_id);
            else
                fprintf(output, "__s->smc_privindex_write(%s, %s, %d, %d, %s, %s, %s, %s, NULL, -1, \"%s\", %d)", input_index, str_string(arrayname), tree->size, tree->sizeexp, output_result, str_string(arraysize1), str_string(arraysize2), buf1, type, tree->thread_id);
        }
    }

    free(type);
    str_free(arrayname);
    str_free(arraysize1);
    str_free(arraysize2);
}

void ast_print_private_indexed_stmt(astexpr tree, char *input_index, char *output_result, int batch_index) {
    char *type = (char *)malloc(sizeof(char) * buffer_size);

    str arrayname = Str("");
    str arraysize1 = Str("");
    str arraysize2 = Str("");

    if (tree->ftype == 0)
        sprintf(type, "int");
    else if (tree->ftype == 1)
        sprintf(type, "float");
    // one-dimension
    if (tree->left->type == IDENT) {
        ast_expr_print(arrayname, tree->left);
        ast_expr_print(arraysize1, tree->left->arraysize);
        indent();
        // for batch private indexing
        if (batch_index != -1)
            fprintf(output, "__s->smc_privindex_read(%s, %s, %s, %s, _picco_batch_counter%d, \"%s\", %d);\n", input_index, str_string(arrayname), output_result, str_string(arraysize1), batch_index, type, tree->thread_id);
        // for singular private indexing
        else
            fprintf(output, "__s->smc_privindex_read(%s, %s, %s, %s, \"%s\", %d);\n", input_index, str_string(arrayname), output_result, str_string(arraysize1), type, tree->thread_id);
    }
    // two-dimension
    else {
        if (tree->left->type == ARRAYIDX) {
            ast_expr_print(arrayname, tree->left->left);
            ast_expr_print(arraysize1, tree->left->arraysize);
            ast_expr_print(arraysize2, tree->left->left->arraysize);
            indent();
            if (batch_index != -1)
                fprintf(output, "__s->smc_privindex_read(%s, %s, %s, %s, %s, _picco_batch_counter%d, \"%s\", %d);\n", input_index, str_string(arrayname), output_result, str_string(arraysize1), str_string(arraysize2), batch_index, type, tree->thread_id);
            else
                fprintf(output, "__s->smc_privindex_read(%s, %s, %s, %s, %s, \"%s\", %d);\n", input_index, str_string(arrayname), output_result, str_string(arraysize1), str_string(arraysize2), type, tree->thread_id);
        }
    }

    free(type);
    str_free(arrayname);
    str_free(arraysize1);
    str_free(arraysize2);
}

void ast_batch_print_stmt_operator(str op, int batch_index, int *private_index, astexpr tree) {

    if (!is_private_indexed(tree)) {
        if (tree->left->type == IDENT) {
            ast_expr_print(op, tree->left);
        } else if (tree->left->type == ARRAYIDX) {
            while (tree->type != IDENT)
                tree = tree->left;
            ast_expr_print(op, tree);
        }
    } else {
        (*private_index)++;
        char *name1 = (char *)malloc(sizeof(char) * buffer_size);
        char *name2 = (char *)malloc(sizeof(char) * buffer_size);

        sprintf(name1, "_picco_private_indexed_inputarray%d", *private_index);
        sprintf(name2, "_picco_private_indexed_outputarray%d", *private_index);

        astexpr e0 = String(name2);
        ast_expr_print(op, e0);

        ast_print_private_indexed_stmt(tree, name1, name2, batch_index);

        free(name1);
        free(name2);
    }
}

void ast_batch_print_stmt_dimension(str dim, astexpr tree) {
    // if the tree is privately indexed, its dimension will be set to 0
    if (is_private_indexed(tree)) {
        str_truncate(dim);
        char *dimension = "0";
        astexpr e = String(dimension);
        ast_expr_print(dim, e);
    } else {
        // two-dimensional op
        // A[i][j]
        if (tree->left->type == ARRAYIDX) {
            str_truncate(dim);
            ast_expr_print(dim, tree->left->left->arraysize);
        }
        // A[i] where A is a two-dimensional array
        else if (tree->arraysize != NULL) {
            str_truncate(dim);
            ast_expr_print(dim, tree->arraysize);
        }
    }
}

void ast_batch_print_stmt_BOP(astexpr tree, int batch_index, int *narray_element_index, int *private_index, char *op, str leftop, str rightop, str leftdim, str rightdim) {
    if ((tree->left->type == ARRAYIDX && tree->left->flag == PRI) && (tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        ast_batch_print_stmt_operator(leftop, batch_index, private_index, tree->left);
        ast_batch_print_stmt_operator(rightop, batch_index, private_index, tree->right);
        ast_batch_print_stmt_dimension(leftdim, tree->left);
        ast_batch_print_stmt_dimension(rightdim, tree->right);
        strcpy(op, BOP_symbols[tree->opid]);
    } else if ((tree->left->type == ARRAYIDX && tree->left->flag == PRI) && !(tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        (*narray_element_index)++;
        str_printf(rightop, "_picco_batch_tmp_array%d", *narray_element_index);
        ast_batch_print_stmt_operator(leftop, batch_index, private_index, tree->left);
        ast_batch_print_stmt_dimension(leftdim, tree->left);
        strcpy(op, BOP_symbols[tree->opid]);
    } else if (!(tree->left->type == ARRAYIDX && tree->left->flag == PRI) && (tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        (*narray_element_index)++;
        str_printf(leftop, "_picco_batch_tmp_array%d", *narray_element_index);
        ast_batch_print_stmt_operator(rightop, batch_index, private_index, tree->right);
        ast_batch_print_stmt_dimension(rightdim, tree->right);
        strcpy(op, BOP_symbols[tree->opid]);
    } else if (!(tree->left->type == ARRAYIDX && tree->left->flag == PRI) && !(tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        (*narray_element_index)++;
        str_printf(leftop, "_picco_batch_tmp_array%d", *narray_element_index);
        strcpy(op, ASS_symbols[0]);
    }
}

int ast_compute_expression_type(astexpr tree) {
    if (tree->type == ASS && tree->right->type == BOP) {
        if (tree->left->ftype == 1 || tree->right->left->ftype == 1 || tree->right->right->ftype == 1)
            return 1;
        else
            return 0;
    } else if (tree->type == BOP) {
        if (tree->left->ftype == 1 || tree->right->ftype == 1)
            return 1;
        else
            return 0;
    } else
        return tree->ftype;
}
void ast_compute_var_size(astexpr tree, char *var_size) {
    int ftype = ast_compute_expression_type(tree);
    int alen_sig = -1, alen_exp = -1, blen_sig = -1, blen_exp = -1;
    int resultlen_sig = -1, resultlen_exp = -1;

    if ((tree->type == ASS && tree->opid == ASS_eq)) {
        astexpr left = tree->right->left;
        astexpr right = tree->right->right;
        if (tree->type == ASS && tree->right->type == BOP) {
            if ((left->type == ARRAYIDX && left->flag == PRI) && (right->type == ARRAYIDX && right->flag == PRI)) {
                alen_sig = left->size;
                alen_exp = left->sizeexp;
                blen_sig = right->size;
                blen_exp = right->sizeexp;
            } else if ((left->type == ARRAYIDX && left->flag == PRI) && (!(right->type == ARRAYIDX) || right->flag == PUB)) {
                alen_sig = left->size;
                alen_exp = left->sizeexp;
                if (right->flag == PRI) {
                    blen_sig = right->size;
                    blen_exp = right->sizeexp;
                }
            } else if ((right->type == ARRAYIDX && right->flag == PRI) && (!(left->type == ARRAYIDX) || left->flag == PUB)) {
                blen_sig = right->size;
                blen_exp = right->sizeexp;
                if (left->flag == PRI) {
                    alen_sig = left->size;
                    alen_exp = left->sizeexp;
                }
            } else {
                if (left->flag == PRI && right->flag == PUB) {
                    alen_sig = left->size;
                    alen_exp = left->sizeexp;
                } else if (left->flag == PUB && right->flag == PRI) {
                    alen_sig = right->size;
                    alen_exp = right->sizeexp;
                } else if (left->flag == PRI && right->flag == PRI) {
                    alen_sig = 32;
                    alen_exp = 9;
                } else {
                    alen_sig = -1;
                    alen_exp = -1;
                }
            }
        } else {
            alen_sig = tree->right->size;
            alen_exp = tree->right->sizeexp;
        }
    }

    else if ((tree->type == ASS && tree->opid != ASS_eq) || tree->type == BOP) {
        alen_sig = tree->left->size;
        alen_exp = tree->left->sizeexp;
        blen_sig = tree->right->size;
        blen_exp = tree->right->sizeexp;
        if (tree->type == BOP)
            if (tree->left->ftype == 1 || tree->right->ftype == 1)
                ftype = 1;
    }
    if (tree->type == ASS) {
        resultlen_sig = tree->left->size;
        resultlen_exp = tree->left->sizeexp;
    }

    if (ftype == 0)
        sprintf(var_size, "%d, %d, %d", alen_sig, blen_sig, resultlen_sig);
    else
        sprintf(var_size, "%d, %d, %d, %d, %d, %d", alen_sig, alen_exp, blen_sig, blen_exp, resultlen_sig, resultlen_exp);
}

void ast_batch_print_stmt(aststmt tree, int batch_index, int statement_index, int *narray_element_index, int *private_index, branchnode current) {
    str leftop, rightop, assignop;
    str leftdim, rightdim, assigndim;
    char *buf1 = (char *)malloc(sizeof(char) * buffer_size);
    char *buf2 = (char *)malloc(sizeof(char) * buffer_size);
    char *buf3 = (char *)malloc(sizeof(char) * buffer_size);
    char *name = (char *)malloc(sizeof(char) * buffer_size);
    char *op = (char *)malloc(sizeof(char) * buffer_size);
    char *type = (char *)malloc(sizeof(char) * buffer_size);
    char *var_size = (char *)malloc(sizeof(char) * buffer_size);

    int privindex = 0;

    leftop = Str("");
    rightop = Str("");
    assignop = Str("");

    leftdim = Str("0");
    rightdim = Str("0");
    assigndim = Str("0");

    int length = 0;
    int private_selection_index = 0;
    int count = 0;
    if (control_sequence_stack_length(private_selection_stack) != 0) {
        private_selection_index = private_selection_stack->head->index;
        count = private_selection_stack->head->batch_index;
    }
    // the batch stmt is not surrounded by any private conditions
    if (private_selection_index == 0 && if_branchnode_height(current) == 0) {
        sprintf(buf1, "NULL");
        sprintf(buf2, "NULL");
        sprintf(buf3, "-1");
    }
    // the batch statement is surrounded by non-batch private conditions
    else if (private_selection_index == 0 && if_branchnode_height(current) != 0) {
        sprintf(buf1, "_picco_condtmp%d", if_branchnode_height(current));
        sprintf(buf2, "NULL");
        sprintf(buf3, "-1");
    }
    // the batch statement is surrounded by batch private conditions
    else {
        sprintf(buf1, "NULL");
        sprintf(buf2, "_picco_batch_array%d", private_selection_index);
        sprintf(buf3, "_picco_batch_counter%d", count);
    }
    // compute the variable size of array elements and store it in var_size

    ast_compute_var_size(tree->u.expr, var_size);
    if (ast_compute_expression_type(tree->u.expr))
        sprintf(type, "float");
    else
        sprintf(type, "int");

    if (tree->u.expr->type == ASS && tree->u.expr->opid == ASS_eq) {
        // for smc_set
        if (is_private_indexed(tree->u.expr->left)) {
            (*private_index)++;
            privindex = (*private_index);
            sprintf(name, "_picco_private_indexed_outputarray%d", privindex);
            astexpr e = String(name);
            ast_expr_print(assignop, e);
        } else {
            ast_batch_print_stmt_operator(assignop, batch_index, private_index, tree->u.expr->left);
            ast_batch_print_stmt_dimension(assigndim, tree->u.expr->left);
        }

        if (tree->u.expr->right->flag == PUB || (tree->u.expr->right->type == CASTEXPR && tree->u.expr->right->left->flag == PUB)) {
            (*narray_element_index)++;
            str_printf(leftop, "_picco_batch_tmp_array%d", *narray_element_index);
            // previous version had " %s ", which would insert whitespace around only the equality operator
            // When op is written inside the smc_batch string below, it caused caused problems  in smc-compute since " = " was unrecognized
            // replaced with "%s" solves the issue
            sprintf(op, "%s", ASS_symbols[tree->u.expr->opid]); 
        } else if (tree->u.expr->right->flag == PRI) {
            if (tree->u.expr->right->type != BOP) {
                // consider the casting
                astexpr tree1 = tree->u.expr->right->type == CASTEXPR ? tree->u.expr->right->left : tree->u.expr->right;
                if (tree1->type == ARRAYIDX) {
                    ast_batch_print_stmt_operator(leftop, batch_index, private_index, tree1);
                    ast_batch_print_stmt_dimension(leftdim, tree1);
                    sprintf(op, " %s ", ASS_symbols[tree->u.expr->opid]);
                } else if (tree1->type != ARRAYIDX) {
                    (*narray_element_index)++;
                    str_printf(leftop, "_picco_batch_tmp_array%d", *narray_element_index);
                    sprintf(op, " %s ", ASS_symbols[tree->u.expr->opid]);
                }
            } else if (tree->u.expr->right->type == BOP)
                ast_batch_print_stmt_BOP(tree->u.expr->right, batch_index, narray_element_index, private_index, op, leftop, rightop, leftdim, rightdim);
        }

    }
    // for ass_add, ass_sub, ass_mul, ass_div, ass_mod
    else if (tree->u.expr->type == ASS && tree->u.expr->opid != ASS_eq) {
        astexpr e = BinaryOperator(BOP_add, tree->u.expr->left, tree->u.expr->right);
        ast_batch_print_stmt_BOP(e, batch_index, narray_element_index, private_index, op, leftop, rightop, leftdim, rightdim);
        strncpy(op, ASS_symbols[tree->u.expr->opid], strlen(ASS_symbols[tree->u.expr->opid]) - 1);
        op[strlen(ASS_symbols[tree->u.expr->opid]) - 1] = '\0';

        if (tree->u.expr->left->type == ARRAYIDX) {
            if (is_private_indexed(tree->u.expr->left)) {
                (*private_index)++;
                privindex = (*private_index);
                sprintf(name, "_picco_private_indexed_outputarray%d", privindex);
                astexpr e = String(name);
                ast_expr_print(assignop, e);
            } else {
                ast_batch_print_stmt_operator(assignop, batch_index, private_index, tree->u.expr->left);
                ast_batch_print_stmt_dimension(assigndim, tree->u.expr->left);
            }
        }
    }
    // for if condition result
    if (tree->u.expr->type == BOP)
        ast_batch_print_private_condition(tree, narray_element_index, private_index, batch_index, statement_index, length, op, leftop, rightop, leftdim, rightdim, current);
    else {
        indent();
        if (!strcmp(str_string(rightop), "") && strcmp(str_string(leftop), "")) {
            if (tree->u.expr->right->type != CASTEXPR)
                fprintf(output, "__s->smc_batch(%s, %s, %s, %s, %s, -1, %s, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, \"%s\", \"%s\", %d);\n", str_string(leftop), str_string(leftop), str_string(assignop), var_size, str_string(leftdim), str_string(assigndim), buf1, buf2, buf3, statement_index, batch_index, op, type, tree->u.expr->thread_id);
            // deal with casting
            else {
                /* INT2FL*/
                if (tree->u.expr->left->ftype == 1 && tree->u.expr->right->left->ftype == 0) {
                    if (tree->u.expr->left->flag == PRI)
                        fprintf(output, "__s->smc_batch_int2fl(%s, %s, %s, %s, %d, %d, %d, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, %d);\n", str_string(leftop), str_string(assignop), str_string(leftdim), str_string(assigndim), tree->u.expr->right->left->size, tree->u.expr->left->size, tree->u.expr->left->sizeexp, buf1, buf2, buf3, statement_index, batch_index, tree->u.expr->thread_id);

                }
                /* FL2INT */
                else if (tree->u.expr->left->ftype == 0 && tree->u.expr->right->left->ftype == 1) {
                    if (tree->u.expr->left->flag == PRI)
                        fprintf(output, "__s->smc_batch_fl2int(%s, %s, %s, %s, %d, %d, %d, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, %d);\n", str_string(leftop), str_string(assignop), str_string(leftdim), str_string(assigndim), tree->u.expr->right->left->size, tree->u.expr->right->left->sizeexp, tree->u.expr->left->size, buf1, buf2, buf3, statement_index, batch_index, tree->u.expr->thread_id);
                }
                /* INT2INT*/
                else if (tree->u.expr->left->ftype == 0 && tree->u.expr->right->left->ftype == 0) {
                    if (tree->u.expr->left->flag == PRI)
                        fprintf(output, "__s->smc_batch_int2int(%s, %s, %s, %s, %d, %d, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, %d);\n", str_string(leftop), str_string(assignop), str_string(leftdim), str_string(assigndim), tree->u.expr->right->left->size, tree->u.expr->left->size, buf1, buf2, buf3, statement_index, batch_index, tree->u.expr->thread_id);

                }
                /* FL2FL */
                else if (tree->u.expr->left->ftype == 1 && tree->u.expr->right->left->ftype == 1) {
                    if (tree->u.expr->left->flag == PRI)
                        fprintf(output, "__s->smc_batch_fl2fl(%s, %s, %s, %s, %d, %d, %d, %d, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, %d);\n", str_string(leftop), str_string(assignop), str_string(leftdim), str_string(assigndim), tree->u.expr->right->left->size, tree->u.expr->right->left->sizeexp, tree->u.expr->left->size, tree->u.expr->left->sizeexp, buf1, buf2, buf3, statement_index, batch_index, tree->u.expr->thread_id);
                }
            }
        }
        // for other BOP
        else if (strcmp(str_string(rightop), "") && strcmp(str_string(leftop), ""))
            fprintf(output, "__s->smc_batch(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, \"%s\", \"%s\", %d);\n", str_string(leftop), str_string(rightop), str_string(assignop), var_size, str_string(leftdim), str_string(rightdim), str_string(assigndim), buf1, buf2, buf3, statement_index, batch_index, op, type, tree->u.expr->thread_id);
        // for BOP with no array element as an operator
        else if (!strcmp(str_string(rightop), "") && !strcmp(str_string(leftop), ""))
            fprintf(output, "__s->smc_batch(NULL, NULL, %s, %s, %s, %s, %s, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, \"%s\", %d);\n", str_string(assignop), var_size, str_string(leftdim), str_string(rightdim), str_string(assigndim), buf1, buf2, buf3, statement_index, batch_index, type, tree->u.expr->thread_id);
        // write back the privately indexed assignment op
        if (is_private_indexed(tree->u.expr->left)) {
            char *input_index = (char *)malloc(sizeof(char) * buffer_size);
            char *output_result = (char *)malloc(sizeof(char) * buffer_size);
            sprintf(input_index, "_picco_private_indexed_inputarray%d", privindex);
            sprintf(output_result, "_picco_private_indexed_outputarray%d", privindex);
            ast_print_private_indexed_stmt_write(tree->u.expr->left, input_index, output_result, batch_index, tree->u.expr->opid, buf1, buf2, buf3);
            fprintf(output, ";\n");
            free(input_index);
            free(output_result);
        }
    }
    free(buf1);
    free(buf2);
    free(buf3);
    free(name);
    free(op);
    free(type);
    str_free(leftop);
    str_free(rightop);
    str_free(assignop);
    str_free(leftdim);
    str_free(rightdim);
    str_free(assigndim);
}

void ast_batch_print_private_condition(aststmt tree, int *narray_element_index, int *private_index, int batch_index, int statement_index, int length, char *op, str leftop, str rightop, str leftdim, str rightdim, branchnode current) {
    char *type = (char *)malloc(sizeof(char) * buffer_size);
    char *buf1 = (char *)malloc(sizeof(char) * buffer_size);
    char *buf2 = (char *)malloc(sizeof(char) * buffer_size);
    char *buf3 = (char *)malloc(sizeof(char) * buffer_size);
    char *leftop_name = (char *)malloc(sizeof(char) * buffer_size);
    char *rightop_name = (char *)malloc(sizeof(char) * buffer_size);
    char *var_size = (char *)malloc(sizeof(char) * buffer_size);

    if (ast_compute_expression_type(tree->u.expr))
        sprintf(type, "float");
    else
        sprintf(type, "int");

    str_truncate(leftop);
    str_truncate(rightop);

    int selection_index = 0;
    int count = 0;
    if (control_sequence_stack_length(private_selection_stack) > 1) {
        selection_index = private_selection_stack->head->next->index;
        count = private_selection_stack->head->next->batch_index;
    }
    // the batch stmt is not surrounded by any private conditions
    if ((if_branchnode_height(current) == control_sequence_stack_length(private_selection_stack)) && control_sequence_stack_length(private_selection_stack) == 1) {
        sprintf(buf1, "NULL");
        sprintf(buf2, "NULL");
        sprintf(buf3, "-1");
    }
    // the batch statement is surrounded by non-batch private conditions
    else if (if_branchnode_height(current) > control_sequence_stack_length(private_selection_stack) && control_sequence_stack_length(private_selection_stack) == 1) {
        sprintf(buf1, "_picco_condtmp%d", if_branchnode_height(current) - control_sequence_stack_length(private_selection_stack));
        sprintf(buf2, "NULL");
        sprintf(buf3, "-1");
    }
    // the batch statement is surrounded by batch private conditions
    else if (control_sequence_stack_length(private_selection_stack) > 1) {
        sprintf(buf1, "NULL");
        sprintf(buf2, "_picco_batch_array%d", selection_index);
        sprintf(buf3, "_picco_batch_counter%d", count);
    }

    int private_selection_index = private_selection_stack->head->index;

    ast_compute_var_size(tree->u.expr, var_size);

    if (current->key == 0)
        ast_batch_print_stmt_BOP(tree->u.expr, batch_index, narray_element_index, private_index, op, leftop, rightop, leftdim, rightdim);
    indent();
    if (!strcmp(str_string(leftop), ""))
        sprintf(leftop_name, "_picco_private_indexed_outputarray%d", *private_index);
    else
        sprintf(leftop_name, " %s ", str_string(leftop));

    if (!strcmp(str_string(rightop), ""))
        sprintf(rightop_name, "_picco_private_indexed_outputarray%d", *private_index);
    else
        sprintf(rightop_name, " %s ", str_string(rightop));

    if (current->key == 0)
        fprintf(output, "__s->smc_batch(%s, %s, _picco_batch_array%d, %s, %s, %s, 0, %s, %s, %s, _picco_batch_index_array%d, _picco_batch_counter%d, \"%s\", \"%s\", %d); \n", leftop_name, rightop_name, private_selection_index, var_size, str_string(leftdim), str_string(rightdim), buf1, buf2, buf3, statement_index, batch_index, op, type, tree->u.expr->thread_id);
    else {
        if (control_sequence_stack_length(private_selection_stack) == 1)
            fprintf(output, "__s->smc_batch(1, _picco_batch_array%d, _picco_batch_array%d, %s, %s, %s, NULL, _picco_batch_counter%d, \"-\", %d);\n", current->if_index, private_selection_index, buf1, buf2, buf3, batch_index, tree->u.expr->thread_id);
        else
            fprintf(output, "__s->smc_batch(_picco_batch_array%d, _picco_batch_array%d, _picco_batch_array%d, %s, %s, %s, NULL, _picco_batch_counter%d, \"-\", %d);\n", current->parent->private_index, current->if_index, private_selection_index, buf1, buf2, buf3, batch_index, tree->u.expr->thread_id);
    }

    free(type);
    free(buf1);
    free(buf2);
    free(buf3);
    free(leftop_name);
    free(rightop_name);
}

// thread_id does not matter in computing private array index as it is local computation
void ast_print_private_index(astexpr tree, char *name) {
    if (tree->left->type == IDENT) {
        astexpr e0 = String(name);
        e0->ftype = tree->right->ftype;
        e0->flag = tree->right->flag;
        e0->isptr = tree->right->isptr;
        astexpr e1 = Assignment(e0, ASS_eq, tree->right);
        e1->ftype = e0->ftype;
        e1->thread_id = tree->thread_id;
        indent();
        ast_expr_show(e1);
        fprintf(output, ";\n");
    } else if (tree->left->type == ARRAYIDX) {
        if (tree->left->left->type == IDENT) {
            // the first or both dimensions are privately indexed
            if (tree->left->right->flag == PRI) {
                /*
        astexpr e0 = tree->left->left->arraysize;
                astexpr e1 = BinaryOperator(BOP_mul, tree->left->right, e0);
                e1->flag = PRI;
                e1->index = 1;
                indent();
                ast_expr_show(e1);
        */
                if (tree->left->right->index > 0)
                    ast_expr_show(tree->left->right);
                str str1 = Str("");
                str str2 = Str("");
                ast_expr_print(str1, tree->left->left->arraysize);
                ast_expr_print(str2, tree->left->right);
                indent();
                indent();
                if (tree->left->right->index > 0)
                    fprintf(output, "__s->smc_mult(_picco_tmp%d, %s, %s, -1, -1, -1, \"int\", %d);\n", tree->left->right->index, str_string(str1), name, tree->thread_id);
                else
                    fprintf(output, "__s->smc_mult(%s, %s, %s, -1, -1, -1, \"int\", %d);\n", str_string(str2), str_string(str1), name, tree->thread_id);

                str_free(str1);
                str_free(str2);
                indent();
                if (tree->right->flag == PRI) {
                    if (tree->right->index > 0) {
                        ast_expr_show(tree->right);
                        indent();
                        indent();
                        fprintf(output, "__s->smc_add(_picco_tmp%d, %s, %s, -1, -1, -1, \"int\", %d);\n", tree->right->index, name, name, tree->thread_id);
                    } else {
                        indent();
                        indent();
                        fprintf(output, "__s->smc_add(");
                        ast_expr_show(tree->right);
                        fprintf(output, ", %s, %s, -1, -1, -1, \"int\", %d);\n", name, name, tree->thread_id);
                    }

                } else if (tree->right->flag == PUB) {
                    indent();
                    fprintf(output, "__s->smc_add(%s, ", name);
                    ast_expr_show(tree->right);
                    fprintf(output, ", %s, -1, -1, -1, \"int\", %d);\n", name, tree->thread_id);
                }
            }
            // the second dimension is privately indexed
            else if (tree->left->right->flag == PUB) {
                astexpr e0 = tree->left->left->arraysize;
                astexpr e1 = UnaryOperator(UOP_paren, tree->left->right);
                astexpr e2 = BinaryOperator(BOP_mul, e0, e1);
                e0->flag = PUB;
                e1->flag = PUB;
                e2->flag = PUB;
                indent();
                indent();
                fprintf(output, "ss_set_ui(%s, ", name);
                ast_expr_show(e2);
                fprintf(output, ");\n");
                indent();
                if (tree->right->index > 0) {
                    ast_expr_show(tree->right);
                    indent();
                    indent();
                    fprintf(output, "__s->smc_add(_picco_tmp%d, %s, %s, -1, -1, -1, \"int\", %d);\n", tree->right->index, name, name, tree->thread_id);
                } else {
                    indent();
                    indent();
                    fprintf(output, "__s->smc_add(");
                    ast_expr_show(tree->right);
                    fprintf(output, ", %s, %s, -1, -1, -1, \"int\", %d);\n", name, name, tree->thread_id);
                }
            }
        }
    }
}

void ast_batch_print_index_operator(str op, astexpr tree, int batch_index, int *private_index) {
    // one-dimension
    if (!is_private_indexed(tree)) {
        if (tree->left->type == IDENT) {
            ast_expr_print(op, tree->right);
        } else if (tree->left->type == ARRAYIDX) {
            // two-dimension
            if (tree->left->left->type == IDENT) {
                astexpr e0 = UnaryOperator(UOP_paren, tree->left->left->arraysize);
                astexpr e1 = UnaryOperator(UOP_paren, tree->left->right);
                astexpr e2 = BinaryOperator(BOP_mul, e0, e1);
                astexpr e3 = BinaryOperator(BOP_add, e2, tree->right);
                ast_expr_print(op, e3);
            }
        }
    } else {
        (*private_index)++;
        char *name = (char *)malloc(sizeof(char) * buffer_size);
        sprintf(name, "_picco_private_indexed_inputarray%d[_picco_ind%d]", *private_index, batch_index);
        ast_print_private_index(tree, name);
        free(name);
    }
}

void ast_batch_print_index_BOP(astexpr tree, int *narray_element_index, int statement_index, int batch_index, str leftop, str rightop, int *private_index) {
    if ((tree->left->type == ARRAYIDX && tree->left->flag == PRI) && (tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        ast_batch_print_index_operator(leftop, tree->left, batch_index, private_index);
        ast_batch_print_index_operator(rightop, tree->right, batch_index, private_index);
    } else if ((tree->left->type == ARRAYIDX && tree->left->flag == PRI) && (!(tree->right->type == ARRAYIDX) || tree->right->flag == PUB)) {
        delete_tmp_array[3 * (*narray_element_index)] = tree->right->flag;
        delete_tmp_array[3 * (*narray_element_index) + 1] = tree->right->ftype;
        delete_tmp_array[3 * (*narray_element_index) + 2] = batch_index;
        (*narray_element_index)++;
        char *name = (char *)malloc(sizeof(char) * buffer_size);
        sprintf(name, "_picco_batch_tmp_array%d[_picco_ind%d]", *narray_element_index, batch_index);
        astexpr e0 = String(name);
        e0->ftype = tree->right->ftype;
        e0->flag = tree->right->flag;
        e0->size = tree->left->size;
        e0->sizeexp = tree->left->sizeexp;
        astexpr e1 = Assignment(e0, ASS_eq, tree->right);
        e1->thread_id = tree->thread_id;
        e1->ftype = e0->ftype;
        indent();
        ast_expr_show(e1);
        fprintf(output, ";\n");
        ast_batch_print_index_operator(leftop, tree->left, batch_index, private_index);
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d+1] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
        free(name);
    }

    else if ((!(tree->left->type == ARRAYIDX) || tree->left->flag == PUB) && (tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        delete_tmp_array[3 * (*narray_element_index)] = tree->left->flag;
        delete_tmp_array[3 * (*narray_element_index) + 1] = tree->left->ftype;
        delete_tmp_array[3 * (*narray_element_index) + 2] = batch_index;
        (*narray_element_index)++;
        char *name = (char *)malloc(sizeof(char) * buffer_size);
        sprintf(name, "_picco_batch_tmp_array%d[_picco_ind%d]", *narray_element_index, batch_index);
        astexpr e0 = String(name);
        e0->ftype = tree->left->ftype;
        e0->flag = tree->left->flag;
        e0->size = tree->right->size;
        e0->sizeexp = tree->right->sizeexp;
        astexpr e1 = Assignment(e0, ASS_eq, tree->left);
        e1->thread_id = tree->thread_id;
        e1->ftype = e0->ftype;
        indent();
        ast_expr_show(e1);
        fprintf(output, ";\n");
        ast_batch_print_index_operator(rightop, tree->right, batch_index, private_index);
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
        free(name);
    } else if ((!(tree->left->type == ARRAYIDX) || tree->left->flag == PUB) && (!(tree->right->type == ARRAYIDX) || tree->right->flag == PUB)) {
        if (tree->left->flag == PRI || tree->right->flag == PRI)
            delete_tmp_array[3 * (*narray_element_index)] = PRI;
        else
            delete_tmp_array[3 * (*narray_element_index)] = PUB;
        delete_tmp_array[3 * (*narray_element_index) + 1] = tree->ftype;
        delete_tmp_array[3 * (*narray_element_index) + 2] = batch_index;
        (*narray_element_index)++;
        char *name = (char *)malloc(sizeof(char) * buffer_size);
        sprintf(name, "_picco_batch_tmp_array%d[_picco_ind%d]", *narray_element_index, batch_index);
        astexpr e0 = String(name);
        e0->ftype = tree->ftype;
        e0->flag = tree->flag;
        if (tree->left->flag == PUB && tree->right->flag == PUB) {
            e0->size = -1;
            e0->sizeexp = -1;
        } else if (tree->left->flag == PRI && tree->right->flag == PUB) {
            e0->size = tree->left->size;
            e0->sizeexp = tree->left->sizeexp;
        } else if (tree->left->flag == PUB && tree->right->flag == PRI) {
            e0->size = tree->right->size;
            e0->sizeexp = tree->right->sizeexp;
        } else {
            e0->size = 32;
            e0->sizeexp = 9;
        }
        astexpr e1 = Assignment(e0, ASS_eq, tree);
        e1->thread_id = tree->thread_id;
        e1->ftype = e0->ftype;
        indent();
        ast_expr_show(e1);
        fprintf(output, ";\n");
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
        free(name);
    }
}

void ast_batch_print_index(aststmt tree, int batch_index, int statement_index, int *narray_element_index, int *delete_tmp_array, int *private_index) {
    // only deals with assignment expression
    str leftop, rightop, assignop, tmp;
    leftop = Str("");
    rightop = Str("");
    assignop = Str("");
    tmp = Str("");

    if (tree->u.expr->type == ASS) {
        ast_batch_print_index_operator(assignop, tree->u.expr->left, batch_index, private_index);
        // right operator is public
        if (tree->u.expr->right->flag == PUB || (tree->u.expr->right->type == CASTEXPR && tree->u.expr->right->left->flag == PUB)) {
            astexpr tree1 = tree->u.expr->right->flag == PUB ? tree->u.expr->right : tree->u.expr->right->left;
            delete_tmp_array[3 * (*narray_element_index)] = PUB;
            delete_tmp_array[3 * (*narray_element_index) + 1] = tree1->ftype;
            delete_tmp_array[3 * (*narray_element_index) + 2] = batch_index;
            (*narray_element_index)++;
            indent();
            fprintf(output, "_picco_batch_tmp_array%d[_picco_ind%d] = ", *narray_element_index, batch_index);
            ast_expr_show(tree1);
            fprintf(output, ";\n");
            indent();
            fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
        } else if (tree->u.expr->right->flag == PRI) {
            if (tree->u.expr->right->type != BOP && tree->u.expr->opid == ASS_eq) {
                astexpr tree1 = tree->u.expr->right->type == CASTEXPR ? tree->u.expr->right->left : tree->u.expr->right;
                if (tree1->type == ARRAYIDX)
                    ast_batch_print_index_operator(leftop, tree1, batch_index, private_index);
                else if (tree1->type != ARRAYIDX) {
                    delete_tmp_array[3 * (*narray_element_index)] = tree1->flag;
                    delete_tmp_array[3 * (*narray_element_index) + 1] = tree1->ftype;
                    delete_tmp_array[3 * (*narray_element_index) + 2] = batch_index;
                    (*narray_element_index)++;
                    char *name = (char *)malloc(sizeof(char) * buffer_size);
                    sprintf(name, "_picco_batch_tmp_array%d[_picco_ind%d]", *narray_element_index, batch_index);
                    astexpr e0 = String(name);
                    e0->ftype = tree1->ftype;
                    e0->flag = tree1->flag;
                    e0->size = tree1->size;
                    e0->sizeexp = tree1->sizeexp;
                    astexpr e1 = Assignment(e0, ASS_eq, tree1);
                    e1->thread_id = tree->u.expr->thread_id;
                    indent();
                    ast_expr_show(e1);
                    fprintf(output, ";\n");
                    indent();
                    fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
                    free(name);
                }
            } else if (tree->u.expr->opid != ASS_eq) {
                astexpr e = BinaryOperator(BOP_add, tree->u.expr->left, tree->u.expr->right);
                ast_batch_print_index_BOP(e, narray_element_index, statement_index, batch_index, leftop, rightop, private_index);
            } else if (tree->u.expr->right->type == BOP)
                ast_batch_print_index_BOP(tree->u.expr->right, narray_element_index, statement_index, batch_index, leftop, rightop, private_index);
        }
    }
    // for if condition (it only works for two operator comparison with each of them being an identifier)
    else if (tree->u.expr->type == BOP)
        ast_batch_print_index_BOP(tree->u.expr, narray_element_index, statement_index, batch_index, leftop, rightop, private_index);
    // print the array index
    if (strcmp(str_string(leftop), "")) {
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d] = %s;\n", statement_index, batch_index, str_string(leftop));
    } else if (!strcmp(str_string(leftop), "")) {
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
    }

    if (strcmp(str_string(rightop), "")) {
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d+1] = %s;\n", statement_index, batch_index, str_string(rightop));
    } else if (!strcmp(str_string(rightop), "")) {
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d+1] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
    }
    if (strcmp(str_string(assignop), "")) {
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d+2] = %s;\n", statement_index, batch_index, str_string(assignop));
    } else if (!strcmp(str_string(assignop), "")) {
        indent();
        fprintf(output, "_picco_batch_index_array%d[3*_picco_ind%d+2] = _picco_ind%d;\n", statement_index, batch_index, batch_index);
    }
    str_free(leftop);
    str_free(rightop);
    str_free(assignop);
    str_free(tmp);
}

void ast_batch_print_prefix_index(aststmt tree, int *batch_index, int *statement_index, int *narray_element_index, int *delete_tmp_array, int *private_index) {
    (*batch_index)++;
    control_sequence_push(*batch_index, batch_stack);
    fprintf(output, "{\n");
    indlev++;
    ast_batch_compute_index(tree, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
    indent();
    fprintf(output, "_picco_ind%d++;\n", batch_stack->head->index);
    indlev--;
    indent();
    fprintf(output, "}\n");
    indlev--;
    control_sequence_pop(batch_stack);
}

void ast_batch_compute_index(aststmt tree, int *batch_index, int *statement_index, int *narray_element_index, int *delete_tmp_array, int *private_index) {
    switch (tree->type) {
    case COMPOUND:
        ast_batch_compute_index(tree->body, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
        break;
    case STATEMENTLIST:
        ast_batch_compute_index(tree->u.next, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
        ast_batch_compute_index(tree->body, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
        break;
    case BATCH:
        ast_batch_print_cond(tree);
        ast_batch_print_prefix_index(tree->body, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
        break;
    case SELECTION: {
        if (tree->u.selection.cond->flag == PUB) {
            indlev++;
            indent();
            fprintf(output, "if (");
            ast_expr_show(tree->u.selection.cond);
            fprintf(output, ")\n");
            indent();
            ast_batch_print_prefix_index(tree->body, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
            if (tree->u.selection.elsebody) {
                indlev++;
                indent();
                fprintf(output, "else\n");
                indent();
                ast_batch_print_prefix_index(tree->u.selection.elsebody, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
            }
        }
        if (tree->u.selection.cond->flag == PRI) {
            (*statement_index)++;
            ast_batch_print_index(Expression(tree->u.selection.cond), batch_stack->head->index, *statement_index, narray_element_index, delete_tmp_array, private_index);
            ast_batch_compute_index(tree->body, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);

            if (tree->u.selection.elsebody) {
                (*statement_index)++;
                ast_batch_compute_index(tree->u.selection.elsebody, batch_index, statement_index, narray_element_index, delete_tmp_array, private_index);
            }
        }
    }
    case EXPRESSION:
        switch (tree->u.expr->type) {
        case ASS:
            (*statement_index)++;
            ast_batch_print_index(tree, batch_stack->head->index, *statement_index, narray_element_index, delete_tmp_array, private_index);
            break;
        }
        break;
    default:
        break;
    }
}

// allocate the memory for each batch statement
void ast_batch_allocate_counter() {
    batch_statement tmp = bss->head;
    int size = 0;
    while (tmp != NULL) {
        indent();
        fprintf(output, "int* _picco_batch_index_array%d = (int*)malloc(sizeof(int) * 3 * _picco_batch_counter%d);\n", tmp->statement_index, tmp->batch_index);
        if (tmp->flag == 1)
            size++;
        tmp = tmp->next;
    }

    tmp = bss->head;
    while (tmp != NULL) {
        if (tmp->flag == 1) {
            indent();
            fprintf(output, "priv_int* _picco_batch_array%d = (priv_int*)malloc(sizeof(priv_int) * _picco_batch_counter%d);\n", size, tmp->batch_index);
            indent();
            fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", tmp->batch_index);
            indent();
            fprintf(output, "{\n");
            indlev++;
            indent();
            fprintf(output, "ss_init(_picco_batch_array%d[_picco_i]);\n", size);
            indlev--;
            indent();
            fprintf(output, "}\n");
            size--;
        }
        tmp = tmp->next;
    }
}

void ast_batch_clear_counter(int private_selection_index, int narray_element_index, int *delete_tmp_array) {
    fprintf(output, "\n");
    batch_statement tmp = bss->head;
    int size = 0, index = 0;
    while (tmp != NULL) {
        indent();
        fprintf(output, "free(_picco_batch_index_array%d);\n", tmp->statement_index);
        if (tmp->flag == 1)
            size++;
        tmp = tmp->next;
    }
    tmp = bss->head;
    while (tmp != NULL) {
        if (tmp->flag == 1) {
            indent();
            fprintf(output, "for (int _picco_j = 0; _picco_j < _picco_batch_counter%d; _picco_j++)\n", tmp->batch_index);
            indent();
            fprintf(output, "{\n");
            indlev++;
            indent();
            fprintf(output, "ss_clear(_picco_batch_array%d[_picco_j]);\n", size);
            indlev--;
            indent();
            fprintf(output, "}\n");
            size--;
        }
        tmp = tmp->next;
    }
    for (index = 0; index < narray_element_index; index++) {
        indent();
        if (delete_tmp_array[3 * index] == PUB)
            fprintf(output, "free(_picco_batch_tmp_array%d);\n", index + 1);
        else if (delete_tmp_array[3 * index] == PRI) {
            if (delete_tmp_array[3 * index + 1] == 0) {
                fprintf(output, "for (int _picco_j = 0; _picco_j < _picco_batch_counter%d; _picco_j++)\n", delete_tmp_array[3 * index + 2]);
                indlev++;
                indent();
                fprintf(output, "ss_clear(_picco_batch_tmp_array%d[_picco_j]);\n", index + 1);
                indlev--;
                indent();
                fprintf(output, "free(_picco_batch_tmp_array%d);\n", index + 1);
            } else if (delete_tmp_array[3 * index + 1] == 1) {
                fprintf(output, "for (int _picco_j = 0; _picco_j < _picco_batch_counter%d; _picco_j++){\n", delete_tmp_array[3 * index + 2]);
                indlev++;
                indent();
                fprintf(output, "for (int _picco_k = 0; _picco_k < 4; _picco_k++)\n");
                indlev++;
                indent();
                fprintf(output, "ss_clear(_picco_batch_tmp_array%d[_picco_j][_picco_k]);\n", index + 1);
                indlev--;
                indent();
                fprintf(output, "free(_picco_batch_tmp_array%d[_picco_j]);\n", index + 1);
                indlev--;
                indent();
                fprintf(output, "}\n");
                indent();
                fprintf(output, "free(_picco_batch_tmp_array%d);\n", index + 1);
            }
        }
    }

    // clear the tmp arrays used for private indexing
    batch_private_index node = bpis->head;
    while (node != NULL) {
        ast_batch_free_array_for_private_indexed_element(node->batch_index, node->private_index, node->ftype);
        node = node->next;
    }
    batch_private_index_stack_free(bpis);
}

void ast_batch_declare_counter(aststmt tree, char *name, int batch_index, branchnode current) {
    // declare the counter variables.
    char *buf = (char *)malloc(sizeof(char) * buffer_size);
    for (ind = 1; ind <= batch_index; ind++) {
        sprintf(buf, "%s%d", name, ind);
        astspec qualifier = Declspec(SPEC_public, 0);
        astspec specifier = Declspec(SPEC_int, 0);
        astspec speclist = Speclist_right(qualifier, specifier);
        // declaration initializer
        astdecl decl = IdentifierDecl(Symbol(buf));
        astexpr expr = Constant(strdup("0"));
        astdecl initdecl = InitDecl(decl, expr);
        // combine the specifier and initializer
        aststmt s = Declaration(speclist, initdecl);
        ast_stmt_show(s, current);
        ast_stmt_free(s);
        indent();
    }
    free(buf);
}
void ast_batch_print_counter(aststmt tree, branchnode current, int *batch_index) {
    char *buf = (char *)malloc(sizeof(char) * buffer_size);
    (*batch_index)++;
    fprintf(output, "{\n");
    indlev++;
    indent();
    sprintf(buf, "_picco_batch_counter%d", *batch_index);
    astexpr e = PreOperator(Identifier(Symbol(buf)), UOP_inc);
    e->flag = PUB;
    aststmt s = Expression(e);
    ast_stmt_show(s, current);
    ast_batch_compute_counter(tree, batch_index, current);
    ast_stmt_free(s);
    indlev--;
    indent();
    fprintf(output, "}\n");
    indlev--;
    free(buf);
}

void ast_batch_compute_counter(aststmt tree, int *batch_index, branchnode current) {
    switch (tree->type) {
    case COMPOUND:
        ast_batch_compute_counter(tree->body, batch_index, current);
        break;
    case STATEMENTLIST:
        ast_batch_compute_counter(tree->u.next, batch_index, current);
        ast_batch_compute_counter(tree->body, batch_index, current);
        break;
    case BATCH:
        // print the batch condition
        ast_batch_print_cond(tree);
        ast_batch_print_counter(tree->body, current, batch_index);
        break;
    case SELECTION:
        if (tree->u.selection.cond->flag == PUB) {
            indlev++;
            indent();
            fprintf(output, "if (");
            ast_expr_show(tree->u.selection.cond);
            fprintf(output, ")\n");
            indent();
            ast_batch_print_counter(tree->body, current, batch_index);
            if (tree->u.selection.elsebody) {
                indlev++;
                indent();
                fprintf(output, "else\n");
                indent();
                ast_batch_print_counter(tree->u.selection.elsebody, current, batch_index);
            }
        } else if (tree->u.selection.cond->flag == PRI) {
            ast_batch_compute_counter(tree->body, batch_index, current);
            if (tree->u.selection.elsebody)
                ast_batch_compute_counter(tree->u.selection.elsebody, batch_index, current);
        }
        break;
    default:
        break;
    }
}

void ast_return_type_of_rop(char **type, astexpr right) {
    /*determine if ptr = &int_var or &float_var or &struct_var*/
    if (right->type == UOP && right->opid == UOP_addr) {
        astexpr tree = right->left;
        while (tree->type == UOP && tree->opid == UOP_paren)
            tree = tree->left;
        while (tree->type == ARRAYIDX)
            tree = tree->left;
        if (tree->type == IDENT) {
            if (!tree->isptr) {
                if (tree->ftype == 0)
                    sprintf(*type, "_int");
                else if (tree->ftype == 1)
                    sprintf(*type, "_float");
                else
                    sprintf(*type, "_struct");
            } else
                sprintf(*type, "\0");
        }
    } else
        sprintf(*type, "\0");
}

void ast_ptr_assignment_show(astexpr left, astexpr right, str rightop) {
    char *type = (char *)malloc(sizeof(char) * buffer_size);
    ast_return_type_of_rop(&type, right);
    fprintf(output, "__s->smc_set%s_ptr(", type);
    ast_expr_show(left);
    fprintf(output, ", %s, ", str_string(rightop));

    if (left->u.sym->type == 2)
        fprintf(output, "\"struct\", %d)", left->thread_id);
    else {
        if (left->ftype == 0)
            fprintf(output, "\"int\", %d)", left->thread_id);
        if (left->ftype == 1)
            fprintf(output, "\"float\", %d)", left->thread_id);
    }

    free(type);
}

int is_ptr_assignment(astexpr tree) {
    if (tree->flag == PRI) {
        if (tree->type == UOP && tree->opid == UOP_paren)
            return is_ptr_assignment(tree->left);
        else if (tree->isptr > 0 || (tree->type == UOP && tree->opid == UOP_addr))
            return 1;
        else
            return 0;
    }
    return 0;
}

void strip_off_bracket_and_dereferences(str op, astexpr tree) {
    while (tree->type == UOP && (tree->opid == UOP_paren) || (tree->opid == UOP_star))
        tree = tree->left;
    ast_expr_print(op, tree);
}

int is_ptr_dereferenced(astexpr tree) {
    if (tree->flag == PRI) {
        if (tree->type == UOP && tree->opid == UOP_paren)
            return is_ptr_dereferenced(tree->left);
        else if (tree->type == UOP && tree->opid == UOP_star)
            return 1;
        else
            return 0;
    }
    return 0;
}

astspec ast_get_struct_name(astexpr tree) {
    struct_node node;
    struct_field field;
    if (tree->left->type != PTRFIELD && tree->left->type != DOTFIELD) {
        node = struct_node_lookup(sns, tree->left->u.sym->struct_type->name->name);
        field = struct_field_lookup(node, tree->u.sym->name);
    } else {
        astspec type = ast_get_struct_name(tree->left);
        node = struct_node_lookup(sns, type->name->name);
        field = struct_field_lookup(node, tree->u.sym->name);
    }
    return field->type;
}

int is_private_struct(astexpr tree) {
    return (tree->type == IDENT) && (tree->u.sym->type == 2) && (tree->isptr) && (!struct_node_get_flag(sns, tree->u.sym->struct_type->name->name));
}

int is_private_struct_field(astexpr tree) {
    // peel the square brackets for arrayindex; (e.g., node->a[1][2] ------> node->a)
    while (tree->type == ARRAYIDX || tree->type == DOTFIELD)
        tree = tree->left;
    if (tree->type == PTRFIELD) {
        while (tree->type == PTRFIELD) {
            char *struct_name = (char *)malloc(sizeof(char) * buffer_size);
            if (tree->left->type == PTRFIELD || tree->left->type == DOTFIELD) {
                astspec type = ast_get_struct_name(tree->left);
                sprintf(struct_name, "%s", type->name->name);
            } else {
                sprintf(struct_name, "%s", tree->left->u.sym->struct_type->name->name);
            }
            if (!struct_node_get_flag(sns, struct_name)) {
                return 1;
            }
            tree = tree->left;
            while (tree->type == ARRAYIDX || tree->type == DOTFIELD)
                tree = tree->left;
            free(struct_name);
        }
        return 0;
    } else
        return 0;
}

int is_private_indexed(astexpr tree) {
    if (tree->flag == PRI) {
        // for one-dimension
        if (tree->type == ARRAYIDX && tree->left->type == IDENT && tree->right->flag == PRI)
            return 1;
        // for two-dimension
        else if (tree->type == ARRAYIDX && tree->left->type == ARRAYIDX &&
                 (tree->right->flag == PRI || tree->left->right->flag == PRI))
            return 1;
        else
            return 0;
    } else
        return 0;
}

void ast_batch_declare_array_for_narrayelement_BOP(astexpr tree, int batch_index, int *narray_element_index, int *private_index) {

    if ((tree->left->type == ARRAYIDX && tree->left->flag == PRI) && (tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        if (is_private_indexed(tree->left)) {
            (*private_index)++;
            ast_batch_allocate_array_for_private_indexed_element(tree->left, batch_index, *private_index);
        }

        if (is_private_indexed(tree->right)) {
            (*private_index)++;
            ast_batch_allocate_array_for_private_indexed_element(tree->right, batch_index, *private_index);
        }

    } else if ((tree->left->type == ARRAYIDX && tree->left->flag == PRI) && !(tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        (*narray_element_index)++;
        ast_batch_allocate_array_for_narrayelement(tree->right, batch_index, *narray_element_index);

        if (is_private_indexed(tree->left)) {
            (*private_index)++;
            ast_batch_allocate_array_for_private_indexed_element(tree->left, batch_index, *private_index);
        }
    } else if (!(tree->left->type == ARRAYIDX && tree->left->flag == PRI) && (tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        (*narray_element_index)++;
        ast_batch_allocate_array_for_narrayelement(tree->left, batch_index, *narray_element_index);

        if (is_private_indexed(tree->right)) {
            (*private_index)++;
            ast_batch_allocate_array_for_private_indexed_element(tree->right, batch_index, *private_index);
        }
    }
    // consider them as a whole.
    else if (!(tree->left->type == ARRAYIDX && tree->left->flag == PRI) && !(tree->right->type == ARRAYIDX && tree->right->flag == PRI)) {
        (*narray_element_index)++;
        ast_batch_allocate_array_for_narrayelement(tree, batch_index, *narray_element_index);
    }
}

void ast_batch_declare_array_for_narrayelement(aststmt tree, int *narray_element_index, int *private_index, int *batch_index) {
    // get the statement information
    switch (tree->type) {
    case COMPOUND:
        ast_batch_declare_array_for_narrayelement(tree->body, narray_element_index, private_index, batch_index);
        break;
    case STATEMENTLIST:
        ast_batch_declare_array_for_narrayelement(tree->u.next, narray_element_index, private_index, batch_index);
        ast_batch_declare_array_for_narrayelement(tree->body, narray_element_index, private_index, batch_index);
        break;
    case BATCH:
        // print the batch condition
        (*batch_index)++;
        control_sequence_push(*batch_index, batch_stack);
        ast_batch_declare_array_for_narrayelement(tree->body, narray_element_index, private_index, batch_index);
        control_sequence_pop(batch_stack);
        break;
    case SELECTION:
        if (tree->u.selection.cond->flag == PUB) {
            (*batch_index)++;
            control_sequence_push(*batch_index, batch_stack);
            ast_batch_declare_array_for_narrayelement(tree->body, narray_element_index, private_index, batch_index);
            control_sequence_pop(batch_stack);
            if (tree->u.selection.elsebody) {
                (*batch_index)++;
                control_sequence_push(*batch_index, batch_stack);
                ast_batch_declare_array_for_narrayelement(tree->u.selection.elsebody, narray_element_index, private_index, batch_index);
                control_sequence_pop(batch_stack);
            }
        } else if (tree->u.selection.cond->flag == PRI) {
            ast_batch_declare_array_for_narrayelement_BOP(tree->u.selection.cond, batch_stack->head->index, narray_element_index, private_index);
            ast_batch_declare_array_for_narrayelement(tree->body, narray_element_index, private_index, batch_index);
            if (tree->u.selection.elsebody) {
                ast_batch_declare_array_for_narrayelement(tree->u.selection.elsebody, narray_element_index, private_index, batch_index);
            }
        }
        break;
    case ITERATION:
        ast_batch_declare_array_for_narrayelement(tree->body, narray_element_index, private_index, batch_index);
        break;
    case EXPRESSION:
        if (tree->u.expr->type == ASS) {
            if (is_private_indexed(tree->u.expr->left)) {
                (*private_index)++;
                ast_batch_allocate_array_for_private_indexed_element(tree->u.expr->left, batch_stack->head->index, *private_index);
            }
            // right operator is public
            if (tree->u.expr->right->flag == PUB || (tree->u.expr->right->type == CASTEXPR && tree->u.expr->right->left->flag == PUB)) {
                (*narray_element_index)++;
                astexpr tree1 = tree->u.expr->right->flag == PUB ? tree->u.expr->right : tree->u.expr->right->left;
                ast_batch_allocate_array_for_narrayelement(tree1, batch_stack->head->index, *narray_element_index);
            }
            // right operator is private
            else if (tree->u.expr->right->flag == PRI) {
                if (tree->u.expr->right->type != BOP) {
                    // consider the casting
                    astexpr tree1 = tree->u.expr->right->type == CASTEXPR ? tree->u.expr->right->left : tree->u.expr->right;
                    if (tree1->type != ARRAYIDX) {
                        (*narray_element_index)++;
                        ast_batch_allocate_array_for_narrayelement(tree1, batch_stack->head->index, *narray_element_index);
                    } else {
                        if (tree->u.expr->opid == ASS_eq) {
                            if (is_private_indexed(tree1)) {
                                (*private_index)++;
                                ast_batch_allocate_array_for_private_indexed_element(tree1, batch_stack->head->index, *private_index);
                            }
                        } else {
                            astexpr e = BinaryOperator(BOP_add, tree->u.expr->left, tree->u.expr->right);
                            ast_batch_declare_array_for_narrayelement_BOP(e, batch_stack->head->index, narray_element_index, private_index);
                        }
                    }
                } else
                    ast_batch_declare_array_for_narrayelement_BOP(tree->u.expr->right, batch_stack->head->index, narray_element_index, private_index);
            }
        }
        break;
    default:
        break;
    }
}

void ast_batch_allocate_array_for_private_indexed_element(astexpr tree, int batch_index, int private_index) {
    if (tree->ftype == 0) {
        indent();
        fprintf(output, "priv_int* _picco_private_indexed_inputarray%d = (priv_int*)malloc(sizeof(priv_int) * _picco_batch_counter%d);\n", private_index, batch_index);
        indent();
        fprintf(output, "priv_int* _picco_private_indexed_outputarray%d = (priv_int*)malloc(sizeof(priv_int) * _picco_batch_counter%d);\n", private_index, batch_index);
        indent();
        fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", batch_index);
        indent();
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "ss_init(_picco_private_indexed_inputarray%d[_picco_i]);\n", private_index);
        indent();
        fprintf(output, "ss_init(_picco_private_indexed_outputarray%d[_picco_i]);\n", private_index);
        indlev--;
        indent();
        fprintf(output, "}\n");
        batch_private_index_push(batch_index, private_index, tree->ftype, bpis);
    } else if (tree->ftype == 1) {
        indent();
        fprintf(output, "priv_int* _picco_private_indexed_inputarray%d = (priv_int*)malloc(sizeof(priv_int) * _picco_batch_counter%d);\n", private_index, batch_index);
        indent();
        fprintf(output, "priv_int** _picco_private_indexed_outputarray%d = (priv_int**)malloc(sizeof(priv_int*) * _picco_batch_counter%d);\n", private_index, batch_index);
        indent();
        fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", batch_index);
        indent();
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "ss_init(_picco_private_indexed_inputarray%d[_picco_i]);\n", private_index);
        indent();
        fprintf(output, "_picco_private_indexed_outputarray%d[_picco_i] = (priv_int*)malloc(sizeof(priv_int) * 4);\n", private_index);
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
        indlev++;
        indent();
        fprintf(output, "ss_init(_picco_private_indexed_outputarray%d[_picco_i][_picco_j]);\n", private_index);
        indlev--;
        indlev--;
        indent();
        fprintf(output, "}\n");
        batch_private_index_push(batch_index, private_index, tree->ftype, bpis);
    }
}

void ast_batch_free_array_for_private_indexed_element(int batch_index, int private_index, int ftype) {

    if (ftype == 0) {
        indent();
        fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", batch_index);
        indent();
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "ss_clear(_picco_private_indexed_inputarray%d[_picco_i]);\n", private_index);
        indent();
        fprintf(output, "ss_clear(_picco_private_indexed_outputarray%d[_picco_i]);\n", private_index);
        indlev--;
        indent();
        fprintf(output, "}\n");
        indent();
        fprintf(output, "free(_picco_private_indexed_inputarray%d);\n", private_index);
        indent();
        fprintf(output, "free(_picco_private_indexed_outputarray%d);\n", private_index);
    } else {
        indent();
        fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", batch_index);
        indent();
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "ss_clear(_picco_private_indexed_inputarray%d[_picco_i]);\n", private_index);
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
        indlev++;
        indent();
        fprintf(output, "ss_clear(_picco_private_indexed_outputarray%d[_picco_i][_picco_j]);\n", private_index);
        indlev--;
        indent();
        fprintf(output, "free(_picco_private_indexed_outputarray%d[_picco_i]);\n", private_index);
        indlev--;
        indent();
        fprintf(output, "}\n");
        indent();
        fprintf(output, "free(_picco_private_indexed_inputarray%d);\n", private_index);
        indent();
        fprintf(output, "free(_picco_private_indexed_outputarray%d);\n", private_index);
    }
}

void ast_batch_allocate_array_for_narrayelement(astexpr tree, int batch_index, int narray_element_index) {
    if (tree->flag == PUB) {
        indent();
        if (tree->ftype == 0)
            fprintf(output, "int* _picco_batch_tmp_array%d = (int*)malloc(sizeof(int) * _picco_batch_counter%d);\n", narray_element_index, batch_index);
        else if (tree->ftype == 1)
            fprintf(output, "float* _picco_batch_tmp_array%d = (float*)malloc(sizeof(float) * _picco_batch_counter%d);\n", narray_element_index, batch_index);
    } else if (tree->flag == PRI) {
        indent();
        if (tree->ftype == 0) {
            fprintf(output, "priv_int* _picco_batch_tmp_array%d = (priv_int*)malloc(sizeof(priv_int) * _picco_batch_counter%d);\n", narray_element_index, batch_index);
            indent();
            fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", batch_index);
            indent();
            fprintf(output, "{\n");
            indlev++;
            indent();
            fprintf(output, "ss_init(_picco_batch_tmp_array%d[_picco_i]);\n", narray_element_index);
            indlev--;
            indent();
            fprintf(output, "}\n");
        } else if (tree->ftype == 1) {
            fprintf(output, "priv_int** _picco_batch_tmp_array%d = (priv_int**)malloc(sizeof(priv_int*) * _picco_batch_counter%d);\n", narray_element_index, batch_index);
            indent();
            fprintf(output, "for (int _picco_i = 0; _picco_i < _picco_batch_counter%d; _picco_i++)\n", batch_index);
            indent();
            fprintf(output, "{\n");
            indlev++;
            indent();
            fprintf(output, "_picco_batch_tmp_array%d[_picco_i] = (priv_int*)malloc(sizeof(priv_int) * 4);\n", narray_element_index);
            indent();
            fprintf(output, "for (int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
            indlev++;
            indent();
            fprintf(output, "ss_init(_picco_batch_tmp_array%d[_picco_i][_picco_j]);\n", narray_element_index);
            indlev--;
            indlev--;
            indent();
            fprintf(output, "}\n");
        }
    }
}

void ast_batch_print_cond(aststmt tree) {
    indent();
    fprintf(output, "for (");
    if (tree->u.iteration.init != NULL) {
        if (tree->u.iteration.init->type == EXPRESSION) {
            if (tree->u.iteration.init->u.expr != NULL)
                ast_expr_show(tree->u.iteration.init->u.expr);
        } else /* Declaration */
        {
            ast_spec_show(tree->u.iteration.init->u.declaration.spec);
            if (tree->u.iteration.init->u.declaration.decl) {
                fprintf(output, " ");
                ast_decl_show(tree->u.iteration.init->u.declaration.decl);
            }
        }
        fprintf(output, "; ");
    } else
        fprintf(output, " ; ");
    if (tree->u.iteration.cond != NULL)
        ast_expr_show(tree->u.iteration.cond);
    fprintf(output, "; ");
    if (tree->u.iteration.incr != NULL)
        ast_expr_show(tree->u.iteration.incr);
    fprintf(output, ")\n");
    indlev++;
    indent();
}

void ast_batch_iter_tree(aststmt tree, int *batch_index, int *statement_index, int *private_selection_index) {
    switch (tree->type) {
    case COMPOUND:
        ast_batch_iter_tree(tree->body, batch_index, statement_index, private_selection_index);
        break;
    case STATEMENTLIST:
        ast_batch_iter_tree(tree->u.next, batch_index, statement_index, private_selection_index);
        ast_batch_iter_tree(tree->body, batch_index, statement_index, private_selection_index);
        break;
    case BATCH: {
        // store the batch condition
        (*batch_index)++;
        control_sequence_push(*batch_index, batch_stack);
        iteration iter = create_iteration(tree->u.iteration.init, tree->u.iteration.cond, tree->u.iteration.incr);
        batch_condition_push(iter, bcs);
        ast_batch_iter_tree(tree->body, batch_index, statement_index, private_selection_index);
        control_sequence_pop(batch_stack);
    } break;
    case SELECTION: {
        // for public selection
        if (tree->u.selection.cond->flag == PUB) {
            (*batch_index)++;
            control_sequence_push(*batch_index, batch_stack);
            ast_batch_iter_tree(tree->body, batch_index, statement_index, private_selection_index);
            control_sequence_pop(batch_stack);
            if (tree->u.selection.elsebody) {
                (*batch_index)++;
                control_sequence_push(*batch_index, batch_stack);
                ast_batch_iter_tree(tree->u.selection.elsebody, batch_index, statement_index, private_selection_index);
                control_sequence_pop(batch_stack);
            }
        }
        // for private selection
        if (tree->u.selection.cond->flag == PRI) {
            (*private_selection_index)++;
            (*statement_index)++;
            batch_statement_push(Expression(tree->u.selection.cond), *statement_index, batch_stack->head->index, 1, bss);
            ast_batch_iter_tree(tree->body, batch_index, statement_index, private_selection_index);
            if (tree->u.selection.elsebody) {
                (*private_selection_index)++;
                (*statement_index)++;
                batch_statement_push(Expression(tree->u.selection.cond), *statement_index, batch_stack->head->index, 1, bss);
                ast_batch_iter_tree(tree->u.selection.elsebody, batch_index, statement_index, private_selection_index);
            }
        }
    } break;
    case ITERATION:
        ast_batch_iter_tree(tree->body, batch_index, statement_index, private_selection_index);
        break;
    case EXPRESSION:
        switch (tree->u.expr->type) {
        case ASS:
            (*statement_index)++;
            batch_statement_push(tree, *statement_index, batch_stack->head->index, 0, bss);
            break;
        }
        break;
    default:
        break;
    }
}

// used for one-dimensional array for smc_input and smc_output
void ast_print_smc_low_dim_show(aststmt tree) {
    is_smc_io = 1;
    ast_expr_show(tree->u.smcops.party);
    fprintf(output, ", ");
    if (tree->u.smcops.size1 == NULL) {
        fprintf(output, "&");
        ast_expr_show(tree->u.smcops.variable);
        fprintf(output, ", \"%s\", %d);", SPEC_symbols[tree->u.smcops.type->subtype], tree->u.smcops.variable->thread_id);
    } else {
        ast_expr_show(tree->u.smcops.variable);
        fprintf(output, ", ");
        ast_expr_show(tree->u.smcops.size1);
        fprintf(output, ", \"%s\", %d);", SPEC_symbols[tree->u.smcops.type->subtype], tree->u.smcops.variable->thread_id);
    }
}

// used for high-dimensional array for smc_input and smc_output
void ast_print_smc_high_dim_show(aststmt tree) {
    fprintf(output, "for(int _picco_i = 0; _picco_i < ");
    ast_expr_show(tree->u.smcops.size1);
    fprintf(output, "; _picco_i++)\n");
    indlev++;
    indent();

    if (tree->subtype == SINPUT)
        fprintf(output, "__s->smc_input(");
    else
        fprintf(output, "__s->smc_output(");

    ast_expr_show(tree->u.smcops.party);
    fprintf(output, ", ");
    ast_expr_show(tree->u.smcops.variable);
    fprintf(output, "[_picco_i], ");
    ast_expr_show(tree->u.smcops.size2);
    fprintf(output, ", \"%s\", %d);", SPEC_symbols[tree->u.smcops.type->subtype], tree->u.smcops.variable->thread_id);
    indlev--;
}

void ast_stmt_smc_show(aststmt tree) {
    switch (tree->subtype) {
    case SINPUT:
        if (tree->u.smcops.size2 == NULL) {
            fprintf(output, "__s->smc_input(");
            ast_print_smc_low_dim_show(tree);
        } else
            ast_print_smc_high_dim_show(tree);
        break;
    case SOUTPUT:
        if (tree->u.smcops.size2 == NULL) {
            fprintf(output, "__s->smc_output(");
            ast_print_smc_low_dim_show(tree);
        } else
            ast_print_smc_high_dim_show(tree);
        break;
    default:
        fprintf(stderr, "[ast_stmt_smc_show]: b u g !!\n");
    }
    fprintf(output, "\n");
}

void ast_stmt_iteration_show(aststmt tree, branchnode current) {
    switch (tree->subtype) {
    case SFOR:
        // print the init
        if (tree->u.iteration.init != NULL) {
            if (tree->u.iteration.init->type == EXPRESSION)
                ast_comma_expr_show(tree->u.iteration.init->u.expr);
            else /* Declaration */
                ast_decl_stmt_show(tree->u.iteration.init, current);
        }
        indent();
        fprintf(output, "for ( ;");
        // print the cond
        if (tree->u.iteration.cond != NULL)
            ast_expr_show(tree->u.iteration.cond);
        fprintf(output, "; )\n");
        indent();
        fprintf(output, "{\n");
        indlev++;
        if (tree->body->type != COMPOUND)
            indent();
        ast_stmt_show(tree->body, current);
        if (tree->u.iteration.incr != NULL) {
            indent();
            ast_comma_expr_show(tree->u.iteration.incr);
        }
        indlev--;
        indent();
        fprintf(output, "}\n");
        break;

    case SWHILE:
        fprintf(output, "while (");
        ast_expr_show(tree->u.iteration.cond);
        fprintf(output, ")\n");
        indent();
        fprintf(output, "{\n");
        indlev++;
        if (tree->body->type != COMPOUND)
            indent();
        ast_stmt_show(tree->body, current);
        indlev--;
        indent();
        fprintf(output, "}\n");
        break;

    case SDO:
        fprintf(output, "do\n");
        indlev++;
        indent();
        ast_stmt_show(tree->body, current);
        indlev--;
        indent();
        fprintf(output, "while (");
        ast_expr_show(tree->u.iteration.cond);
        fprintf(output, ");\n\n");
        break;

    default:
        fprintf(stderr, "[ast_stmt_iteration_show]: b u g !!\n");
    }
}

void ast_stmt_selection_show(aststmt tree, branchnode current) {
    switch (tree->subtype) {
    case SSWITCH:
        fprintf(output, "switch (");
        ast_expr_show(tree->u.selection.cond);
        fprintf(output, ")\n");
        indent();
        ast_stmt_show(tree->body, current);
        break;

    case SIF:
        // for pub cond
        if (!ast_check_priv_if(tree->u.selection.cond)) {
            fprintf(output, "if (");
            ast_expr_show(tree->u.selection.cond);
            fprintf(output, ")\n");

            indlev++;
            indent();
            ast_stmt_show(tree->body, current);
            indlev--;

            if (tree->u.selection.elsebody) {
                indent();
                fprintf(output, "else\n");
                indlev++;
                indent();
                ast_stmt_show(tree->u.selection.elsebody, current);
                indlev--;
            }
        }
        // for pri cond
        else if (ast_check_priv_if(tree->u.selection.cond)) {
            astexpr pubcond = NULL;
            astexpr pricond = NULL;
            aststmt ifstmt = NULL;
            branchnode right = NULL;
            branchnode left = NULL;
            priv_if_index++;
            int priv_else_index = priv_if_index;

            if (tree->u.selection.cond->type == COMMALIST) {
                ast_filter_cond(tree->u.selection.cond, &pubcond, &pricond);
                if_push(pricond, if_top);
                ast_expr_show(pricond);
            } else {
                if_push(tree->u.selection.cond, if_top);
                ast_expr_show(tree->u.selection.cond);
            }
            // ms stores the pointers that are assigned values in the body of private if/else
            mvarstack ms = NULL;
            if (tree->u.selection.elsebody) {
                ms = mvar_stack_new();
                ast_iter_tree(tree, ms);
                ast_compute_selection_completeness(tree, ms);
            }
            left = if_branchnode_insert(current, ms, 0, -1, -1, current->current_label, priv_if_index);
            indlev++;
            indent();
            indlev++;
            fprintf(output, "{\n");
            if (pubcond != NULL) {
                ifstmt = If(pubcond, tree->body, NULL);
                ast_if_stmt_show(ifstmt, left, tree->u.selection.cond->thread_id);
            } else
                ast_if_stmt_show(tree->body, left, tree->u.selection.cond->thread_id);
            indlev--;
            indent();
            indlev--;
            fprintf(output, "}\n");

            if (tree->u.selection.elsebody) {
                priv_if_index++;
                indlev++;
                indent();
                indlev++;
                fprintf(output, "{\n");
                right = if_branchnode_insert(current, NULL, 1, -1, -1, current->current_label, priv_if_index);
                if (tree->u.selection.elsebody->type != COMPOUND)
                    indent();
                ast_if_stmt_show(tree->u.selection.elsebody, right, tree->u.selection.cond->thread_id);
                ast_shrink_ptr(ms, priv_else_index, current->current_label, tree->u.selection.cond->thread_id);
                indlev--;
                indent();
                indlev--;
                fprintf(output, "}\n");
            }
            if_pop(if_top);
            if (left != NULL)
                if_branchtree_remove(current, left);
            if (right != NULL)
                if_branchtree_remove(current, right);
        }
        break;
    default:
        fprintf(stderr, "[ast_stmt_selection_show]: b u g !!\n");
    }
}

void ast_filter_cond(astexpr cond, astexpr *pubcond, astexpr *pricond) {
    int i = 0;
    int j = 0;
    while (1) {
        if (cond->right->flag == PUB) {
            if (i == 0)
                *pubcond = cond->right;
            else
                *pubcond = Astexpr(COMMALIST, *pubcond, cond->right);
            i++;
        }

        else if (cond->right->flag == PRI) {
            if (j == 0)
                *pricond = cond->right;
            else
                *pricond = BinaryOperator(BOP_land, *pricond, cond->right);
            j++;
        }

        if (cond->left->type != COMMALIST)
            break;
        cond = cond->left;
    }

    if (cond->left->flag == PUB)
        if (i == 0)
            *pubcond = cond->left;
        else
            *pubcond = Astexpr(COMMALIST, *pubcond, cond->left);

    else if (j == 0)
        *pricond = cond->left;
    else
        *pricond = BinaryOperator(BOP_land, *pricond, cond->left);
}

void ast_if_stmt_show(aststmt tree, branchnode current, int thread_id) {
    int line = 0;
    int start_addr = 0;
    int is_priv_int_appear = 0, is_priv_float_appear = 0;

    if (if_branchnode_height(current) == 1) {
        ast_compute_type_of_assignment_var(tree, &is_priv_int_appear, &is_priv_float_appear);
        if (is_priv_int_appear)
            ast_tmp_decl_show("_picco_priv_int_", 1, 1);
        if (is_priv_float_appear)
            ast_float_tmp_decl_show("_picco_priv_float_", 1, 1);
    }
    ast_compute_condition(if_top->head, current, thread_id);
    ast_stmt_show(tree, current);
    ast_tmp_clear_show("_picco_condtmp", if_branchnode_height(current), if_branchnode_height(current));

    // clear
    if (if_branchnode_height(current) == 1) {
        if (is_priv_int_appear)
            ast_tmp_clear_show("_picco_priv_int_tmp", 1, 1);
        if (is_priv_float_appear)
            ast_float_tmp_clear_show("_picco_priv_float_tmp", 1, 1);
    }
}

void ast_compute_type_of_assignment_var(aststmt tree, int *is_priv_int_appear, int *is_priv_float_appear) {
    switch (tree->type) {
    case COMPOUND:
        ast_compute_type_of_assignment_var(tree->body, is_priv_int_appear, is_priv_float_appear);
        break;
    case STATEMENTLIST:
        ast_compute_type_of_assignment_var(tree->u.next, is_priv_int_appear, is_priv_float_appear);
        ast_compute_type_of_assignment_var(tree->body, is_priv_int_appear, is_priv_float_appear);
        break;
    case SELECTION:
        ast_compute_type_of_assignment_var(tree->body, is_priv_int_appear, is_priv_float_appear);
        if (tree->u.selection.elsebody)
            ast_compute_type_of_assignment_var(tree->u.selection.elsebody, is_priv_int_appear, is_priv_float_appear);
        break;
    case ITERATION:
        ast_compute_type_of_assignment_var(tree->body, is_priv_int_appear, is_priv_float_appear);
        break;
    case EXPRESSION:
        if (tree->u.expr->type == ASS) {
            if (tree->u.expr->left->ftype == 0)
                *is_priv_int_appear = 1;
            else
                *is_priv_float_appear = 1;
        }
        break;
    default:
        break;
    }
}

void ast_is_selection_complete(aststmt tree, astexpr var, int *result) {
    str name = Str("");
    ast_expr_print(name, var);
    switch (tree->type) {
    case COMPOUND:
        ast_is_selection_complete(tree->body, var, result);
        break;
    case STATEMENTLIST:
        ast_is_selection_complete(tree->u.next, var, result);
        ast_is_selection_complete(tree->body, var, result);
        break;
    case SELECTION: {
        int left = 0;
        int right = 0;
        ast_is_selection_complete(tree->body, var, &left);
        if (tree->u.selection.elsebody)
            ast_is_selection_complete(tree->u.selection.elsebody, var, &right);
        if (left + right == 2)
            *result = 1;
        break;
    }
    case ITERATION:
        ast_is_selection_complete(tree->body, var, result);
        break;
    case EXPRESSION:
        switch (tree->u.expr->type) {
        case ASS: {
            str tmp_name;
            if (tree->u.expr->left->type == IDENT) {
                if (tree->u.expr->left->isptr > 0) {
                    tmp_name = Str("");
                    ast_expr_print(tmp_name, tree->u.expr->left);
                    if (!strcmp(str_string(name), str_string(tmp_name)))
                        *result = 1;
                    str_free(tmp_name);
                }
            } else if (tree->u.expr->left->type == ARRAYIDX) {
                // modify this if start dealing with array of pointers
            }
            break;
        }
        }
        break;
    default:
        break;
    }
    str_free(name);
}

void ast_shrink_ptr(mvarstack mvartable, int priv_else_index, int parent_index, int thread_id) {
    if (mvartable != NULL) {
        mvar node = mvartable->head;
        while (node != NULL) {
            if (node->is_complete == 1) {
                str name = Str("");
                ast_expr_print(name, node->var_name);
                indent();
                fprintf(output, "__s->smc_shrink_ptr(%s, %d, %d, %d);\n", str_string(name), priv_else_index, parent_index, thread_id);
                node->is_complete = 0;
                str_free(name);
            }
            node = node->next;
        }
    }
}

void ast_compute_selection_completeness(aststmt tree, mvarstack mvartable) {
    mvar node = mvartable->head;
    int left = 0, right = 0;
    while (node != NULL) {
        left = 0;
        right = 0;
        ast_is_selection_complete(tree->body, node->var_name, &left);
        ast_is_selection_complete(tree->u.selection.elsebody, node->var_name, &right);
        if (left == 1 && right == 1)
            node->is_complete = 1;
        node = node->next;
    }
}

void ast_iter_tree(aststmt tree, mvarstack mvartable) {
    switch (tree->type) {
    case COMPOUND:
        ast_iter_tree(tree->body, mvartable);
        break;
    case STATEMENTLIST:
        ast_iter_tree(tree->u.next, mvartable);
        ast_iter_tree(tree->body, mvartable);
        break;
    case SELECTION:
        ast_iter_tree(tree->body, mvartable);
        if (tree->u.selection.elsebody)
            ast_iter_tree(tree->u.selection.elsebody, mvartable);
        break;
    case ITERATION:
        ast_iter_tree(tree->body, mvartable);
        break;
    case EXPRESSION:
        switch (tree->u.expr->type) {
        case ASS:
            if (tree->u.expr->left->type == IDENT) {
                if (tree->u.expr->left->isptr > 0)
                    mvar_push(tree->u.expr->left, mvartable);
            } else if (tree->u.expr->left->type == ARRAYIDX) {
                // modify this if start dealing with array of pointers
            }
            break;
        }
        break;
    default:
        break;
    }
}

void ast_compute_condition(condnode cnode, branchnode bnode, int thread_id) {
    int index = cnode->element->index;
    int height = if_branchnode_height(bnode);
    indent();
    fprintf(output, "priv_int _picco_condtmp%d;\n", height);
    indent();
    fprintf(output, "ss_init(_picco_condtmp%d);\n", height);
    indent();

    if (bnode->key == 0)
        fprintf(output, "ss_set(_picco_condtmp%d, _picco_tmp%d);\n", height, index);
    else
        fprintf(output, "__s->smc_sub(1, _picco_tmp%d, _picco_condtmp%d, 1, 1, 1, \"int\", %d);\n", index, height, thread_id);

    if (height != 1) {
        indent();
        fprintf(output, "__s->smc_mult(_picco_condtmp%d, _picco_condtmp%d, _picco_condtmp%d, 1, 1, 1, \"int\", %d);\n", height - 1, height, height, thread_id);
    }
}

void ast_stmt_labeled_show(aststmt tree, branchnode current) {
    switch (tree->subtype) {
    case SLABEL:
        fprintf(output, "%s :\n", tree->u.label->name);
        break;
    case SCASE:
        fprintf(output, "case ");
        ast_expr_show(tree->u.expr);
        fprintf(output, " :\n");
        break;
    case SDEFAULT:
        fprintf(output, "default :\n");
        break;
    default:
        fprintf(stderr, "[ast_stmt_labeled_show]: b u g !!\n");
    }
    indlev++;
    indent();
    ast_stmt_show(tree->body, current);
    indlev--;
}

void ast_free_memory_for_local_variables(ltablelist tablelist) {
    ltable table = ltablelist_pop(tablelist);
    if (table == NULL)
        return;
    lvar var = table->head;
    while (var != NULL) {
        ast_handle_memory_for_private_variable(var->decl, var->spec, "", 1);
        var = var->next;
    }
    ltable_free(table);
}

/**
 * Display/Traverses abstract syntax tree statements and handle output streams.
 * - Processes each statement node in the tree based on its type.
 * - Handles expression statements, including private variables and batch loops.
 * - Prints compound statements (blocks enclosed in braces).
 * - Manages memory cleanup after compound statements.
 * - Handles declaration statements and temporary variable declaration.
 *
 * @param tree The AST representing statements to display.
 * @param output_filename The output file where statements will be displayed.
 **/

void ast_stmt_show(aststmt tree, branchnode current) {
    if (tree->file != NULL && tree->file != _curfile)
        if (tree->type != FUNCDEF && tree->type != STATEMENTLIST) {
            _curfile = tree->file;
            if (strcmp("injected_code", tree->file->name)) {
                fprintf(output, "# %d \"%s\"\n", tree->l, tree->file->name);
                indent();
            }
        }
    switch (tree->type) {
    case JUMP:
        ast_stmt_jump_show(tree, current);
        break;
    case ITERATION:
        ast_stmt_iteration_show(tree, current);
        break;
    case SELECTION:
        ast_stmt_selection_show(tree, current);
        break;
    case LABELED:
        ast_stmt_labeled_show(tree, current);
        break;
    case EXPRESSION:
        if (tree->u.expr != NULL) {
            // expr is within private-if but outside the batch loop
            if (global_batch_flag != 1 && if_branchnode_height(current) != 0) {
                if (tree->u.expr->type == ASS || tree->u.expr->type == POSTOP || tree->u.expr->type == PREOP)
                    ast_priv_expr_show(tree->u.expr, current, tree->gflag); // send the gflag to this function
            } else {
                ast_expr_show(tree->u.expr);
                fprintf(output, ";\n");
            }
        }
        break;
    case COMPOUND:
        indent();
        fprintf(output, "{\n");
        if (tree->body) {
            ltablelist_push(current->tablelist);
            indlev++;
            indent();
            ast_stmt_show(tree->body, current); /* Ends in \n */
            indlev--;
            if (is_return_void || (if_branchnode_height(current) != 0 || ltablelist_length(current->tablelist) > 1)) {
                if (if_branchnode_height(current) == 0 && ltablelist_length(current->tablelist) <= 1) {
                    ast_clear_all_tmp_variables(current);
                    is_return_void = 0;
                } else {
                    ast_free_memory_for_local_variables(current->tablelist);
                }
            }
        }
        indent();
        fprintf(output, "}\n");
        break;
    case STATEMENTLIST: {
        aststmt ch;
        int lastdef = 0;

        /* If my rightmost child of the left subtree is a DECLARATION and
         * the leftmost child of my right subtree is not, then this is
         * where declarations end.
         */
        for (ch = tree->u.next; ch->type == STATEMENTLIST; ch = ch->body)
            ;
        if (ch->type == DECLARATION) {
            for (ch = tree->body; ch->type == STATEMENTLIST; ch = ch->u.next)
                ;
            if (ch->type != DECLARATION)
                lastdef = 1;
        }
        if (declared == 0 && enterfunc == 1) {
            ast_temporary_variable_declaration();
            declared = 1;
        }
        ast_stmt_show(tree->u.next, current);
        /* 1 empty line after the last declaration */
        fprintf(output, "\n");
        indent();
        ast_stmt_show(tree->body, current);
        break;
    }
    case DECLARATION:
        if (!tree->is_stmt_for_sng)
            ast_decl_stmt_show(tree, current);
        else
            ast_decl_sng_stmt_show(tree);
        break;
    case VERBATIM:
        fprintf(output, "%s\n", tree->u.code);
        break;
    case FUNCDEF:
        fprintf(output, "\n"); /* Make it stand out */
        indent();
        tmp_index = tree->num_tmp;
        tmp_float_index = tree->num_float_tmp;
        is_priv_int_ptr_appear = tree->is_priv_int_ptr_appear;
        is_priv_float_ptr_appear = tree->is_priv_float_ptr_appear;
        is_priv_int_index_appear = tree->is_priv_int_index_appear;
        is_priv_float_index_appear = tree->is_priv_float_index_appear;
        is_priv_int_struct_field_appear = tree->is_priv_int_struct_field_appear;
        is_priv_float_struct_field_appear = tree->is_priv_float_struct_field_appear;
        is_priv_struct_ptr_struct_field_appear = tree->is_priv_struct_ptr_struct_field_appear;
        is_priv_int_ptr_struct_field_appear = tree->is_priv_int_ptr_struct_field_appear;
        is_priv_float_ptr_struct_field_appear = tree->is_priv_float_ptr_struct_field_appear;
        contain_priv_if_flag = tree->contain_priv_if_flag;

        declared = 0;
        enterfunc = 1;

        // check for void return type
        if ((tree->u.declaration.spec->subtype == SPEC_Rlist && tree->u.declaration.spec->u.next->subtype == SPEC_void) || tree->u.declaration.spec->subtype == SPEC_void)
            is_return_void = 1;
        ast_spec_show(tree->u.declaration.spec);
        fprintf(output, " ");
        ast_decl_show(tree->u.declaration.decl);
        fprintf(output, "\n");
        if (tree->u.declaration.dlist) {
            indlev++;
            indent();
            ast_stmt_show(tree->u.declaration.dlist, current);
            indlev--;
        }
        indent();
        ast_stmt_show(tree->body, current);

        is_priv_int_ptr_appear = 0;
        is_priv_float_ptr_appear = 0;
        is_priv_int_index_appear = 0;
        is_priv_float_index_appear = 0;
        is_priv_int_struct_field_appear = 0;
        is_priv_float_struct_field_appear = 0;
        is_priv_struct_ptr_struct_field_appear = 0;
        is_priv_int_ptr_struct_field_appear = 0;
        is_priv_float_ptr_struct_field_appear = 0;
        fprintf(output, "\n");
        enterfunc = 0; /* Make it stand out */
        break;
    case SMC:
        ast_stmt_smc_show(tree);
        break;
    case BATCH:
        global_batch_flag++;
        ast_stmt_batch_show(tree, current);
        global_batch_flag--;
        break;
    case OMPSTMT:
        ast_ompcon_show(tree->u.omp, current);
        fprintf(output, "\n");
        break;
    case OX_STMT:
        ast_oxcon_show(tree->u.ox, current);
        fprintf(output, "\n");
        break;
    default:
        fprintf(stderr, "[ast_stmt_show]: b u g !!\n");
    }
}

void ast_expr_priv_index(astexpr tree) {
    char *input_index = (char *)malloc(sizeof(char) * buffer_size);
    char *output_result = (char *)malloc(sizeof(char) * buffer_size);

    if (is_private_indexed(tree->left)) {
        sprintf(input_index, "_picco_priv_ind1");
        if (tree->ftype == 0)
            sprintf(output_result, "_picco_priv_tmp1");
        else if (tree->ftype == 1)
            sprintf(output_result, "_picco_priv_ftmp1");
        ast_print_private_index(tree->left, input_index);
        ast_print_private_indexed_stmt(tree->left, input_index, output_result, -1);
    }
    if (is_private_indexed(tree->right)) {
        sprintf(input_index, "_picco_priv_ind2");
        if (tree->ftype == 0)
            sprintf(output_result, "_picco_priv_tmp2");
        else if (tree->ftype == 1)
            sprintf(output_result, "_picco_priv_ftmp2");

        ast_print_private_index(tree->right, input_index);
        ast_print_private_indexed_stmt(tree->right, input_index, output_result, -1);
    }

    free(input_index);
    free(output_result);
}

void ast_compute_dereferences_and_pointers(astexpr e, int *dereferences, int *pointers) {
    *dereferences = 0;
    *pointers = 0;

    while (e->type != IDENT) {
        if (e->type == UOP) {
            if (e->opid == UOP_star)
                (*dereferences)++;
            e = e->left;
        }
    }
    *pointers = e->isptr;
}

void ast_print_pi_ptr_operator(astexpr tree, int index) {
    char *output_result;
    int dereferences = 0, pointers = 0;
    if (is_ptr_dereferenced(tree)) {
        ast_compute_dereferences_and_pointers(tree, &dereferences, &pointers);
        if (dereferences == pointers) {
            if (tree->ftype == 1)
                output_result = "_picco_priv_ftmp";
            else if (tree->ftype == 0)
                output_result = "_picco_priv_tmp";
        } else {
            if (tree->ftype == 1)
                output_result = "_picco_tmp_float_ptr";
            else if (tree->ftype == 0)
                output_result = "_picco_tmp_int_ptr";
        }
        fprintf(output, "%s%d, ", output_result, index);
    } else if (is_private_indexed(tree)) {
        if (tree->ftype == 1)
            output_result = "_picco_priv_ftmp";
        else if (tree->ftype == 0)
            output_result = "_picco_priv_tmp";
        fprintf(output, "%s%d, ", output_result, index);
    } else if (is_private_struct_field(tree)) {
        int ispointer, field_type;
        char *var = (char *)malloc(sizeof(char) * buffer_size);
        char *struct_name = (char *)malloc(sizeof(char) * buffer_size);
        astexpr tree1 = NULL;

        if (tree->type == ARRAYIDX && tree->left->type == PTRFIELD)
            tree1 = tree->left;
        else if (tree->type == ARRAYIDX && tree->left->type == ARRAYIDX)
            tree1 = tree->left->left;
        else
            tree1 = tree;
        if (tree1->left->type == PTRFIELD || tree1->left->type == DOTFIELD) {
            astspec type = ast_get_struct_name(tree1->left);
            sprintf(struct_name, "%s", type->name->name);
        } else {
            sprintf(struct_name, "%s", tree1->left->u.sym->struct_type->name->name);
        }
        ast_return_struct_field_info(struct_name, tree1->u.sym->name, &ispointer, &field_type);
        ast_declare_temp_for_struct_field(&var, ispointer, field_type, index);
        fprintf(output, "%s, ", var);
        free(struct_name);
        free(var);
    }
}

void ast_init_temp_for_struct_field(char **var, int ispointer, int field_type) {
    if (ispointer) {
        indent();
        fprintf(output, "__s->smc_clear_ptr(&(%s));\n", *var);
    } else {
        if (field_type == 0) {
            indent();
            fprintf(output, "ss_set_ui(%s, 0);\n", *var);
        } else if (field_type == 1) {
            indent();
            fprintf(output, "for(int _picco_ind = 0; _picco_ind < 4; _picco_ind++)\n");
            indlev++;
            indent();
            fprintf(output, "ss_set_ui(%s[_picco_ind], 0);\n", *var);
            indlev--;
        }
    }
}

void ast_declare_temp_for_struct_field(char **var, int ispointer, int field_type, int id) {
    if (ispointer) {
        if (field_type == 0)
            sprintf(*var, "_picco_str_field_tmp_int_ptr%d", id);
        else if (field_type == 1)
            sprintf(*var, "_picco_str_field_tmp_float_ptr%d", id);
        else if (field_type == 2)
            sprintf(*var, "_picco_str_field_tmp_struct_ptr%d", id);
    } else {
        if (field_type == 0)
            sprintf(*var, "_picco_str_field_tmp_int%d", id);
        else if (field_type == 1)
            sprintf(*var, "_picco_str_field_tmp_float%d", id);
    }
}

void ast_expr_refer_single_struct_field_helper(char *var_name, char *struct_name, char *field_name, int *id, int is_address, int read_or_write, int private_if_index, char *rightop, astexpr tree, int thread_id) {
    char *priv_cond = (char *)malloc(sizeof(char) * buffer_size);
    char *array_index = (char *)malloc(sizeof(char) * buffer_size);

    int ispointer, field_type;
    str first_dimension = Str(""), second_dimension = Str("");
    if (private_if_index == -1)
        sprintf(priv_cond, "NULL");
    else
        sprintf(priv_cond, "_picco_condtmp%d", private_if_index);

    if (tree) {
        if (tree->type == ARRAYIDX && tree->left->type == IDENT) {
            ast_expr_print(first_dimension, tree->right);
            sprintf(array_index, "%s, ", str_string(first_dimension));
        } else if (tree->type == ARRAYIDX && tree->left->type == ARRAYIDX) {
            ast_expr_print(first_dimension, tree->left->right);
            ast_expr_print(second_dimension, tree->right);
            sprintf(array_index, "%s, %s, ", str_string(first_dimension), str_string(second_dimension));
        }
    } else
        sprintf(array_index, "");

    if (!struct_node_get_flag(sns, struct_name)) {
        char *var = (char *)malloc(sizeof(char) * buffer_size);
        char *var_type = (char *)malloc(sizeof(char) * buffer_size);
        ast_return_struct_field_info(struct_name, field_name, &ispointer, &field_type);
        if (field_type == 0)
            sprintf(var_type, "int");
        else if (field_type == 1)
            sprintf(var_type, "float");
        else
            sprintf(var_type, "struct");

        if (!read_or_write) {
            if (*id <= 3) {
                if (field_type == 2)
                    *id += 3;
            } else
                *id -= 3;
            ast_declare_temp_for_struct_field(&var, ispointer, field_type, *id);
            ast_init_temp_for_struct_field(&var, ispointer, field_type);
        }

        // call struct field dereference function
        indent();
        if (!read_or_write)
            fprintf(output, "%s_%s(%s, %s, %s%d, %s, %d);\n", struct_name, field_name, var_name, var, array_index, read_or_write, priv_cond, thread_id);
        else {
            // if the right operator is an address of variable, we create a temporary ptr variable, and store it.
            if (is_address) {
                fprintf(output, "__s->smc_clear_ptr(&(_picco_tmp_%s_ptr));\n", var_type);
                indent();
                fprintf(output, "__s->smc_set_%s_ptr(_picco_tmp_%s_ptr, %s, \"%s\", %d);\n", var_type, var_type, rightop, var_type, thread_id);
                indent();
                fprintf(output, "%s_%s(%s, _picco_tmp_%s_ptr, %s%d, %s, %d);\n", struct_name, field_name, var_name, var_type, array_index, read_or_write, priv_cond, thread_id);
            } else
                fprintf(output, "%s_%s(%s, %s, %s%d, %s, %d);\n", struct_name, field_name, var_name, rightop, array_index, read_or_write, priv_cond, thread_id);
        }
        free(var_type);
        free(var);
    }

    free(priv_cond);
    free(array_index);
    str_free(first_dimension);
    str_free(second_dimension);
}

astspec ast_expr_refer_single_struct_field(astexpr tree, int *id, int is_address, int read_or_write, int private_if_index, char *rightop, int *recursive_level, char **var, int *tag, int thread_id) {
    char *var_name = (char *)malloc(sizeof(char) * buffer_size);
    char *struct_name = (char *)malloc(sizeof(char) * buffer_size);
    char *field_name = (char *)malloc(sizeof(char) * buffer_size);
    struct_node node;
    struct_field field;
    astexpr tree1 = NULL;

    // peel the brackets of array field.
    if (tree->type == ARRAYIDX) {
        tree1 = tree;
        while (tree->type == ARRAYIDX)
            tree = tree->left;
    }

    // multiple field dereference
    if (tree->left != NULL && tree->left->type == PTRFIELD) {
        (*recursive_level)++;
        int struct_id = 0;
        astspec type = ast_expr_refer_single_struct_field(tree->left, id, 0, 0, private_if_index, NULL, recursive_level, var, tag, thread_id);
        (*recursive_level)--;
        sprintf(struct_name, "%s", type->name->name);
        sprintf(field_name, "%s", tree->u.sym->name);

        node = struct_node_lookup(sns, struct_name);
        field = struct_field_lookup(node, field_name);

        if (node->contain_pub_field) {
            strcat(*var, "->");
            strcat(*var, tree->left->u.sym->name);
        } else {
            if (*tag == 1) {
                sprintf(var_name, "%s->%s", *var, tree->left->u.sym->name);
                *tag = 0;
            } else
                sprintf(var_name, "_picco_str_field_tmp_struct_ptr%d", *id);
        }
        // for now we only support static array on non-pointer primitive types
        if (*recursive_level == 0) {
            ast_expr_refer_single_struct_field_helper(var_name, struct_name, field_name, id, is_address, read_or_write, private_if_index, rightop, tree1, thread_id);
            if (*id > 3) {
                *id -= 3;
                char *var1 = (char *)malloc(sizeof(char) * buffer_size);
                char *var2 = (char *)malloc(sizeof(char) * buffer_size);
                int ispointer = 0, field_type = 0;
                ast_return_struct_field_info(struct_name, field_name, &ispointer, &field_type);
                ast_declare_temp_for_struct_field(&var1, ispointer, field_type, *id);
                ast_declare_temp_for_struct_field(&var2, ispointer, field_type, *id + 3);
                if (field_type == 2) {
                    indent();
                    fprintf(output, "__s->smc_clear_ptr(&(%s));\n", var1);
                    indent();
                    fprintf(output, "__s->smc_set_ptr(%s, %s, \"struct\", %d);\n", var1, var2, thread_id);
                }
                free(var1);
                free(var2);
            }
        } else
            ast_expr_refer_single_struct_field_helper(var_name, struct_name, field_name, id, 0, 0, private_if_index, NULL, NULL, thread_id);

        free(var_name);
        free(struct_name);
        free(field_name);
        return field->type;
    }
    // single field dereference
    else if (tree->type == PTRFIELD) {
        // if the variable is a struct type
        if (tree->left->type == IDENT)
            sprintf(struct_name, "%s", tree->left->u.sym->struct_type->name->name);
        else if (tree->left->type == DOTFIELD) {
            astspec s = ast_get_struct_name(tree->left);
            sprintf(struct_name, "%s", s->name->name);
        }
        sprintf(field_name, "%s", tree->u.sym->name);
        str left = Str("");
        node = struct_node_lookup(sns, struct_name);
        field = struct_field_lookup(node, field_name);
        ast_expr_print(left, tree->left);

        if (node->contain_pub_field) {
            *tag = 1;
            sprintf(*var, "%s", str_string(left));
        } else
            sprintf(var_name, "%s", str_string(left));

        ast_expr_refer_single_struct_field_helper(var_name, struct_name, field_name, id, is_address, read_or_write, private_if_index, rightop, tree1, thread_id);
        if (*id > 3) {
            *id -= 3;
            char *var1 = (char *)malloc(sizeof(char) * buffer_size);
            char *var2 = (char *)malloc(sizeof(char) * buffer_size);
            int ispointer = 0, field_type = 0;
            ast_return_struct_field_info(struct_name, field_name, &ispointer, &field_type);
            ast_declare_temp_for_struct_field(&var1, ispointer, field_type, *id);
            ast_declare_temp_for_struct_field(&var2, ispointer, field_type, *id + 3);
            if (field_type == 2) {
                indent();
                fprintf(output, "__s->smc_clear_ptr(&(%s));\n", var1);
                indent();
                fprintf(output, "__s->smc_set_ptr(%s, %s, \"struct\", %d);\n", var1, var2, thread_id);
            }
            free(var1);
            free(var2);
        }

        free(var_name);
        free(struct_name);
        free(field_name);
        str_free(left);
        return field->type;
    }

    free(var_name);
    free(struct_name);
    free(field_name);
}

void ast_expr_refer_struct_field(astexpr tree, int private_if_index) {
    int recursive_level = 0;
    char *var_left = (char *)malloc(sizeof(char) * buffer_size);
    char *var_right = (char *)malloc(sizeof(char) * buffer_size);
    int tag_left = 0, tag_right = 0;
    int id1 = 1, id2 = 2;
    if (is_private_struct_field(tree->left))
        ast_expr_refer_single_struct_field(tree->left, &id1, 0, 0, private_if_index, NULL, &recursive_level, &var_left, &tag_left, tree->left->thread_id);
    recursive_level = 0;
    if (is_private_struct_field(tree->right))
        ast_expr_refer_single_struct_field(tree->right, &id2, 0, 0, private_if_index, NULL, &recursive_level, &var_right, &tag_right, tree->right->thread_id);
    free(var_left);
    free(var_right);
}

void ast_expr_ptr_dereference(astexpr tree, int private_if_index) {
    char *input_index = (char *)malloc(sizeof(char) * buffer_size);
    char *output_result = (char *)malloc(sizeof(char) * buffer_size);
    char *priv_cond = (char *)malloc(sizeof(char) * buffer_size);
    int dereferences = 0, pointers = 0;
    str left = Str("");
    str right = Str("");

    if (private_if_index == -1)
        sprintf(priv_cond, "NULL");
    else
        sprintf(priv_cond, "_picco_condtmp%d", private_if_index);

    if (is_ptr_dereferenced(tree->left)) {
        strip_off_bracket_and_dereferences(left, tree->left);
        ast_compute_dereferences_and_pointers(tree->left, &dereferences, &pointers);
        if (tree->ftype == 0) {
            if (dereferences == pointers)
                sprintf(output_result, "_picco_priv_tmp1");
            else
                sprintf(output_result, "_picco_tmp_int_ptr1");
            indent();
            fprintf(output, "__s->smc_dereference_read_ptr(%s, %s, %d, %s, \"int\", %d);\n", str_string(left), output_result, dereferences, priv_cond, tree->thread_id);
        } else if (tree->ftype == 1) {
            if (dereferences == pointers)
                sprintf(output_result, "_picco_priv_ftmp1");
            else
                sprintf(output_result, "_picco_tmp_float_ptr1");
            indent();
            fprintf(output, "__s->smc_dereference_read_ptr(%s, %s, %d, %s, \"float\", %d);\n", str_string(left), output_result, dereferences, priv_cond, tree->thread_id);
        }
    }

    if (is_ptr_dereferenced(tree->right)) {
        strip_off_bracket_and_dereferences(right, tree->right);
        ast_compute_dereferences_and_pointers(tree->right, &dereferences, &pointers);
        if (tree->ftype == 0) {
            if (dereferences == pointers)
                sprintf(output_result, "_picco_priv_tmp2");
            else
                sprintf(output_result, "_picco_tmp_int_ptr2");
            indent();
            fprintf(output, "__s->smc_dereference_read_ptr(%s, %s, %d, %s, \"int\", %d);\n", str_string(right), output_result, dereferences, priv_cond, tree->thread_id);
        } else if (tree->ftype == 1) {
            if (dereferences == pointers)
                sprintf(output_result, "_picco_priv_ftmp2");
            else
                sprintf(output_result, "_picco_tmp_float_ptr2");
            indent();
            fprintf(output, "__s->smc_dereference_read_ptr(%s, %s, %d, %s, \"float\", %d);\n", str_string(right), output_result, dereferences, priv_cond, tree->thread_id);
        }
    }
    free(input_index);
    free(output_result);
    free(priv_cond);
    str_free(left);
    str_free(right);
}

/*by ghada*/
void ast_expr_pmalloc_show(astexpr tree) {
    /*The implementation depends on the data type of the private data
     *         we have three cases:
     *          case 1: private int
     *          case 2: private float
     *          case 3: struct
     *                                         */
    switch (tree->u.dtype->spec->type) {
    case SPEC:
        if (tree->u.dtype->spec->subtype == SPEC_int) {
            indent();
            fprintf(output, "_picco_temp_ = malloc(");
            ast_expr_show(tree->left);
            fprintf(output, "*sizeof(priv_int));\n");
            indent();
            fprintf(output, "for (int _picco_i = 0; _picco_i < ");
            ast_expr_show(tree->left);
            fprintf(output, "; _picco_i++)\n");
            indlev++;
            indent();
            fprintf(output, "ss_init(_picco_temp_[_picco_i]);\n");
            indlev--;
        } else if (tree->u.dtype->spec->subtype == SPEC_float) {
            indent();
            fprintf(output, "_picco_temp_ = malloc(");
            ast_expr_show(tree->left);
            fprintf(output, "*sizeof(priv_int*));\n");
            indent();
            fprintf(output, "for (int _picco_i = 0; _picco_i < ");
            ast_expr_show(tree->left);
            fprintf(output, "; _picco_i++){\n");
            indlev++;
            indent();
            indent();
            fprintf(output, "_picco_temp_[_picco_i] = (priv_int*)malloc(sizeof(priv_int) * 4);\n");
            indent();
            indent();
            fprintf(output, "for (int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
            indlev++;
            indent();
            indent();
            fprintf(output, "ss_init(_picco_temp_[_picco_i][_picco_j]);\n");
            indlev--;
            indlev--;
            indent();
            fprintf(output, "}\n");
        }
        break;
    case SUE:
        switch (tree->u.dtype->spec->subtype) {
        case SPEC_struct:
        case SPEC_union:
            fprintf(output, "_picco_temp_ = malloc(");
            ast_expr_show(tree->left);
            fprintf(output, "*sizeof(struct %s));\n", tree->u.dtype->spec->name->name);
            indent();
            fprintf(output, "for (int _picco_i = 0; _picco_i < ");
            ast_expr_show(tree->left);
            fprintf(output, "; _picco_i++)\n");
            indlev++;
            indent();
            fprintf(output, "%s_init(&((struct %s*)_picco_temp_)[_picco_i]);\n", tree->u.dtype->spec->name->name, tree->u.dtype->spec->name->name);
            indlev--;
            break;
        default:
            fprintf(stderr, "pmalloc can be used with int, float, struct, union only!!\n");
        }
        break;
    case SPECLIST:
        switch (tree->u.dtype->spec->subtype) {
        case SPEC_Rlist:
            tree->u.dtype->spec = tree->u.dtype->spec->u.next;
            ast_expr_pmalloc_show(tree);
            break;
        case SPEC_Llist:
            break;
        case SPEC_enum:
            break;
        default:
            fprintf(stderr, "pmalloc can be used with int, float, struct, union only!!\n");
        }
        break;
    default:
        fprintf(stderr, "pmalloc can be used with int, float, struct, union only!!\n");
    }
}

void ast_priv_expr_show(astexpr tree, branchnode current, int gflag) {
    // implement privtmp = x;
    char *type = NULL;
    char *inttype = "int";
    char *floattype = "float";
    char *op;
    int is_complete = 0;
    if (tree->left->ftype == 0) {
        type = inttype;
        op = "_picco_priv_int_tmp1";
    } else if (tree->left->ftype == 1) {
        type = floattype;
        op = "_picco_priv_float_tmp1";
    }
    arg_str = Str("");
    ast_expr_print(arg_str, tree->left);
    char *left = str_string(arg_str);

    // pointer assignment
    if (tree->left->isptr > 0 && !is_private_struct_field(tree->left)) {
        str right;
        // deal with the case e.g, ptr1 = **ptr2
        if (is_ptr_dereferenced(tree->right)) {
            indent();
            if (tree->right->ftype == 0)
                fprintf(output, "__s->smc_clear_ptr(&_picco_tmp_int_ptr1);\n");
            else if (tree->right->ftype == 1)
                fprintf(output, "__s->smc_clear_ptr(&_picco_tmp_float_ptr1);\n");
            astexpr e0 = Constant(strdup("0"));
            astexpr e1 = BinaryOperator(BOP_add, tree->right, e0);
            e1->ftype = tree->right->ftype;
            e1->thread_id = tree->right->thread_id;
            ast_expr_ptr_dereference(e1, if_branchnode_height(current));
            if (tree->right->ftype == 0)
                right = Str("_picco_tmp_int_ptr1");
            else
                right = Str("_picco_tmp_float_ptr1");
        }
        // deal with the case e.g., ptr1 = ptr2 or &var1
        if (is_ptr_assignment(tree->right) && !is_private_struct_field(tree->right)) {
            right = Str("");
            ast_expr_print(right, tree->right);
        }
        // deal with the case e.g., ptr1 = str->ptr2
        if (is_private_struct_field(tree->right)) {
            astexpr e0 = Constant(strdup("0"));
            astexpr e1 = BinaryOperator(BOP_add, tree->right, e0);
            e1->ftype = tree->right->ftype;
            e1->thread_id = tree->right->thread_id;
            ast_expr_refer_struct_field(e1, if_branchnode_height(current));
            int ispointer, field_type;
            char *struct_name = (char *)malloc(sizeof(char) * buffer_size);
            char *rightop = (char *)malloc(sizeof(char) * buffer_size);
            right = Str("");
            astexpr tree1 = NULL;
            if (tree->right->type == ARRAYIDX && tree->right->left->type == PTRFIELD)
                tree1 = tree->right->left;
            else if (tree->right->type == ARRAYIDX && tree->right->left->type == ARRAYIDX)
                tree1 = tree->right->left->left;
            else
                tree1 = tree->right;
            if (tree1->left->type == PTRFIELD || tree1->left->type == DOTFIELD) {
                astspec type = ast_get_struct_name(tree1->left);
                sprintf(struct_name, "%s", type->name->name);
            } else
                sprintf(struct_name, "%s", tree1->left->u.sym->struct_type->name->name);

            ast_return_struct_field_info(struct_name, tree1->u.sym->name, &ispointer, &field_type);
            ast_declare_temp_for_struct_field(&rightop, ispointer, field_type, 1);
            astexpr name = String(rightop);
            ast_expr_print(right, name);
            free(struct_name);
            free(rightop);
        }
        indent();
        char *right_type = (char *)malloc(sizeof(char) * buffer_size);
        ast_return_type_of_rop(&right_type, tree->right);
        fprintf(output, "__s->smc_update%s_ptr(%s, %s, _picco_condtmp%d, %d, %d);\n", right_type, left, str_string(right), if_branchnode_height(current), priv_if_index, tree->right->thread_id);
        free(right_type);
        str_free(right);
    }
    // pointer dereference
    else if (is_ptr_dereferenced(tree->left) || is_private_indexed(tree->left) || is_private_struct_field(tree->left)) {
        ast_priv_assignment_show(tree, if_branchnode_height(current));
        fprintf(output, ";\n");
    } else {

        /*privtmp = x*/
        indent();
        if (tree->left->ftype == 0) {
            fprintf(output, "__s->smc_set(%s, %s, %d, %d, \"%s\", %d);\n", left, op, tree->left->size, tree->left->size, type, tree->left->thread_id);
        } else if (tree->left->ftype == 1) {
            fprintf(output, "__s->smc_set(%s, %s, %d, %d, %d, %d, \"%s\", %d);\n", left, op, tree->left->size, tree->left->sizeexp, tree->left->size, tree->left->sizeexp, type, tree->left->thread_id);
        }
        /*evaluate expr*/
        indent();
        ast_expr_show(tree);
        fprintf(output, ";\n");
        /*x-privtmp*/
        indent();
        fprintf(output, "__s->smc_priv_eval(%s, %s, _picco_condtmp%d, %d);\n", left, op, if_branchnode_height(current), tree->left->thread_id);
        // if(tree->left->ftype == 0)
        //       fprintf(output, "__s->smc_sub(%s, %s, %s, %d, %d, %d, \"%s\", %d);\n", left, op, left, tree->left->size, tree->left->size, tree->left->size, type, tree->left->thread_id);
        // else if(tree->left->ftype == 1)
        //       fprintf(output, "__s->smc_sub(%s, %s, %s, %d, %d, %d, %d, %d, %d, \"%s\", %d);\n", left, op, left, tree->left->size, tree->left->sizeexp, tree->left->size, tree->left->sizeexp, tree->left->size, tree->left->sizeexp, type, tree->left->thread_id);

        /* (x-privtmp) * cond */
        // indent();
        // if(tree->left->ftype == 0)
        //      fprintf(output, "__s->smc_mult(%s, _picco_condtmp%d, %s, %d, 1, %d, \"%s\", %d);\n", left, if_branchnode_height(current), left, tree->left->size, tree->left->size, type, tree->left->thread_id);
        // else if(tree->left->ftype == 1)
        //       fprintf(output, "__s->smc_mult(%s, _picco_condtmp%d, %s, %d, %d, 1, %d, %d, \"%s\", %d);\n", left, if_branchnode_height(current), left, tree->left->size, tree->left->sizeexp, tree->left->size, tree->left->sizeexp, type, tree->left->thread_id);

        /* (x-privtmp) * cond + privtmp */
        // indent();
        // if(tree->left->ftype == 0)
        //      fprintf(output, "__s->smc_add(%s, %s, %s, %d, %d, %d, \"%s\", %d);\n", left, op, left, tree->left->size, tree->left->size, tree->left->size, type, tree->left->thread_id);
        // else if(tree->left->ftype == 1)
        //      fprintf(output, "__s->smc_add(%s, %s, %s, %d, %d, %d, %d, %d, %d, \"%s\", %d);\n", left, op, left, tree->left->size, tree->left->sizeexp, tree->left->size, tree->left->sizeexp, tree->left->size, tree->left->sizeexp, type, tree->left->thread_id);
    }
}

/**
 * This function recursively traverses the abstract syntax tree (AST) 
 * representing expressions and prints them out according to certain conditions.
 * It handles various types of expressions including identifiers, constant values,
 * function calls, array indexing, struct fields, pointer fields, casts, operators,
 * assignments, designated initializers, and more.
 * 
 * Note: this function is responsible for printing the name and the value of the variables
 * all over the generated code. Hence, there are different flags used to separate private global 
 * init from the rest. 
 * 
 * - Processes each expression in the AST 'tree' based on its type.
 * - Manages printing of identifiers, constant values, strings, function calls, array indices, etc.
 * - Handles special cases like typecasting, struct fields, pointer fields, and function calls.
 * - Determines whether to print expressions into the global stream or the regular output stream.
 * - Handles initialization of private variables and special SMC function calls.
 * - Manages the indentation and formatting of the output.
 * - Supports typecasting between different data types including int and float.
 *
 * @param tree The AST representing expression statements to display.
 **/
void ast_expr_show(astexpr tree) {
    switch (tree->type) {
    case IDENT:
        if (is_smc_set == 1 && is_smc_io == 0 && is_priv == 1) {
            fprintf(global_stream, "%s", tree->u.sym->name);
            // printf("\nName: %s\n", tree->u.sym->name); // the private name gets printed here
        } else {
            fprintf(output, "%s", tree->u.sym->name);
        }
        break;
    case CONSTVAL:
        if (is_smc_set == 1 && is_smc_io == 0 && is_priv == 1) {
            fprintf(global_stream, "%s", tree->u.str);
            // printf("\nValue: %s\n", tree->u.str); // the const value for both private and public gets printed here
        } else {
            fprintf(output, "%s", tree->u.str);
        }
        break;
    case STRING: 
        fprintf(output, "%s", tree->u.str);
        break;
    case FUNCCALL:
        arg_str = Str("");
        ast_expr_print(arg_str, tree->left);

        if (!strcmp(str_string(arg_str), "smcopen"))
            fprintf(output, "%s", "__s->smc_open");
        else
            ast_expr_show(tree->left);

        fprintf(output, "(");
        if (tree->right)
            ast_expr_show(tree->right);

        if (!strcmp(str_string(arg_str), "smcopen"))
            fprintf(output, ", %d)", tree->thread_id);
        else
            fprintf(output, ")");

        break;
    case ARRAYIDX:
        ast_expr_show(tree->left);
        fprintf(output, "[");
        ast_expr_show(tree->right);
        fprintf(output, "]");
        break;
    case DOTFIELD:
        ast_expr_show(tree->left);
        fprintf(output, ".%s", tree->u.sym->name);
        break;
    case PTRFIELD:
        // retrieve the struct type for tree->left
        // if it does contain any public field, we will use its original name
        // otherwise, we will add suffix at the end.
        // if tree->left->u.sym->struct_type == NULL -------> it is the struct declared by openMP constructs

        if (tree->left->u.sym->struct_type == NULL || struct_node_get_flag(sns, tree->left->u.sym->struct_type->name->name)) {
            ast_expr_show(tree->left);
            fprintf(output, "->%s", tree->u.sym->name);
        }
        break;
    case BRACEDINIT:
        if (tree->left->type != COMMALIST)
            fprintf(output, "{ ");
        else {
            fprintf(output, "{\n");
            indlev += 2;
            indent();
        }
        ast_expr_show(tree->left);
        if (tree->left->type != COMMALIST)
            fprintf(output, " }");
        else {
            fprintf(output, "\n");
            indlev--;
            indent();
            indlev--;
            fprintf(output, "}");
        }
        break;
    // NEEDS MORE WORK
    case CASTEXPR:
        /* private typecasting */
        if (tree->u.dtype->spec->subtype == SPEC_float || tree->u.dtype->spec->subtype == SPEC_int || (tree->u.dtype->spec->subtype == SPEC_Rlist && tree->u.dtype->spec->body->subtype == SPEC_private)) {
            /* conversion to float */
            if (tree->u.dtype->spec->subtype == SPEC_float || tree->u.dtype->spec->subtype == SPEC_Rlist &&
                                                                tree->u.dtype->spec->body->subtype == SPEC_private && tree->u.dtype->spec->u.next->subtype == SPEC_float) {
                /* Int2FL */
                if (tree->left->ftype == 0) {
                    ast_priv_single_expr_show(tree->left);
                    fprintf(output, "__s->smc_int2fl(");
                    ast_priv_cast_helper_show(tree->left);
                    fprintf(output, "_picco_ftmp%d, %d, %d, %d, %d);\n", tree->index, tree->left->size, tree->size, tree->sizeexp, tree->thread_id);
                }
                /* FL2FL */
                else {
                    ast_priv_single_expr_show(tree->left);
                    fprintf(output, "__s->smc_fl2fl(");
                    ast_priv_cast_helper_show(tree->left);
                    fprintf(output, "_picco_ftmp%d, %d, %d, %d, %d, %d);\n", tree->index, tree->left->size, tree->left->sizeexp, tree->size, tree->sizeexp, tree->thread_id);
                }
            }
            /* conversion to int */
            if (tree->u.dtype->spec->subtype == SPEC_int || tree->u.dtype->spec->subtype == SPEC_Rlist &&
                                                                tree->u.dtype->spec->body->subtype == SPEC_private && tree->u.dtype->spec->u.next->subtype == SPEC_int) {
                /* Int2Int */
                if (tree->left->ftype == 0) {
                    ast_priv_single_expr_show(tree->left);
                    fprintf(output, "__s->smc_int2int(");
                    ast_priv_cast_helper_show(tree->left);
                    fprintf(output, "_picco_tmp%d, %d, %d, %d);\n", tree->index, tree->left->size, tree->size, tree->thread_id);
                }
                /* FL2Int */
                else {
                    ast_priv_single_expr_show(tree->left);
                    fprintf(output, "__s->smc_fl2int(");
                    ast_priv_cast_helper_show(tree->left);
                    fprintf(output, "_picco_tmp%d, %d, %d, %d, %d);\n", tree->index, tree->left->size, tree->left->sizeexp, tree->size, tree->thread_id);
                }
            }
        } else {
            fprintf(output, "(");
            ast_decl_show(tree->u.dtype);
            fprintf(output, ") ");
            ast_expr_show(tree->left);
        }
        break;
    case CONDEXPR:
        ast_expr_show(tree->u.cond);
        fprintf(output, " ? ");
        ast_expr_show(tree->left);
        fprintf(output, " : ");
        ast_expr_show(tree->right);
        break;

    case UOP:
        if (tree->flag != PRI)
            fprintf(output, "%s", UOP_symbols[tree->opid]);
        if (tree->opid == UOP_sizeoftype || tree->opid == UOP_sizeof)
            fprintf(output, "(");
        if (tree->opid == UOP_sizeoftype || tree->opid == UOP_typetrick)
            ast_decl_show(tree->u.dtype);
        else if (tree->flag == PRI && tree->opid == UOP_lnot) {
            fprintf(output, "__s->smc_eqeq(");
            ast_expr_show(tree->left);
            fprintf(output, ", 0, _picco_tmp%d, %d, -1, %d, \"int\", %d);\n", tree->index, tree->left->size, tree->size, tree->thread_id);
            indent();
        } else
            ast_expr_show(tree->left);
        if (tree->opid == UOP_paren && tree->flag != PRI || tree->opid == UOP_sizeoftype || tree->opid == UOP_sizeof)
            fprintf(output, ")");
        break;

    case BOP:
        /* for pub operation */
        if (tree->left->flag == PUB && tree->right->flag == PUB) {
            ast_expr_show(tree->left);
            fprintf(output, " %s ", BOP_symbols[tree->opid]);
            ast_expr_show(tree->right);
        }
        /* for pri operation */
        else {
            if (tree->left->index > 0) {
                ast_expr_show(tree->left);
            }
            if (tree->right->index > 0) {
                if (tree->left->index > 0)
                    indent();
                ast_expr_show(tree->right);
            }
            // for correct indentation
            if (tree->left->index > 0 || tree->right->index > 0)
                indent();
            /*check if any of the operands is privately indexed*/
            ast_expr_priv_index(tree);
            /*check if any of the operands is pointer dereferenced*/
            ast_expr_ptr_dereference(tree, -1);
            /*check if any of the operands is a struct field referenced by private pointer*/
            ast_expr_refer_struct_field(tree, -1);
            switch (tree->opid) {
            case BOP_add:
                ast_smc_show("__s->smc_add", tree);
                break;
            case BOP_sub:
                ast_smc_show("__s->smc_sub", tree);
                break;
            case BOP_mul:
                ast_smc_show("__s->smc_mult", tree);
                break;
            case BOP_div:
                ast_smc_show("__s->smc_div", tree);
                break;
            case BOP_lt:
                ast_smc_show("__s->smc_lt", tree);
                break;
            case BOP_gt:
                ast_smc_show("__s->smc_gt", tree);
                break;
            case BOP_leq:
                ast_smc_show("__s->smc_leq", tree);
                break;
            case BOP_geq:
                ast_smc_show("__s->smc_geq", tree);
                break;
            case BOP_eqeq:
                ast_smc_show("__s->smc_eqeq", tree);
                break;
            case BOP_neq:
                ast_smc_show("__s->smc_neq", tree);
                break;
            case BOP_lor:
                ast_smc_show("__s->smc_lor", tree);
                break;
            case BOP_land:
                ast_smc_show("__s->smc_land", tree);
                break;
            case BOP_band:
                ast_smc_show("__s->smc_band", tree);
                break;
            case BOP_xor:
                ast_smc_show("__s->smc_xor", tree);
                break;
            case BOP_bor:
                ast_smc_show("__s->smc_bor", tree);
                break;
            case BOP_dot:
                ast_smc_show("__s->smc_dot", tree);
                break;
            case BOP_shl:
                ast_smc_show("__s->smc_shl", tree);
                break;
            case BOP_shr:
                ast_smc_show("__s->smc_shr", tree);
                break;
            }

            // print two operands (consider if the operator is privately indexed)
            if (is_private_indexed(tree->left) || is_ptr_dereferenced(tree->left) || is_private_struct_field(tree->left))
                ast_print_pi_ptr_operator(tree->left, 1);
            else {
                if (tree->left->index <= 0)
                    ast_expr_show(tree->left);
                else {
                    if (tree->ftype == 1)
                        fprintf(output, "_picco_ftmp%d", tree->left->index);
                    else if (tree->ftype == 0)
                        fprintf(output, "_picco_tmp%d", tree->left->index);
                }
                fprintf(output, ", ");
            }

            if (is_private_indexed(tree->right) || is_ptr_dereferenced(tree->right) || is_private_struct_field(tree->right))
                ast_print_pi_ptr_operator(tree->right, 2);
            else {
                if (tree->right->index <= 0)
                    ast_expr_show(tree->right);
                else {
                    if (tree->ftype == 1)
                        fprintf(output, "_picco_ftmp%d", tree->right->index);
                    else if (tree->ftype == 0)
                        fprintf(output, "_picco_tmp%d", tree->right->index);
                }
                fprintf(output, ", ");
            }
            // if the operation is a dot product
            if (tree->opid == BOP_dot) {
                ast_expr_show(tree->left->arraysize);
                fprintf(output, ", ");
            }

            if (tree->left->arraytype != 1 && tree->right->arraytype != 1) {
                if (tree->ftype == 1)
                    fprintf(output, "_picco_ftmp%d, ", tree->index);
                else if (tree->ftype == 0)
                    fprintf(output, "_picco_tmp%d, ", tree->index);
                if (tree->right->ftype == 1 || tree->left->ftype == 1) {
                    if (tree->opid != BOP_lt && tree->opid != BOP_gt && tree->opid != BOP_leq && tree->opid != BOP_geq && tree->opid != BOP_eqeq && tree->opid != BOP_neq)
                        fprintf(output, "%d, %d, %d, %d, %d, %d, \"float\", %d);\n", tree->left->size, tree->left->sizeexp, tree->right->size, tree->right->sizeexp, tree->size, tree->sizeexp, tree->thread_id);
                    else
                        fprintf(output, "%d, %d, %d, %d, 1, \"float\", %d);\n", tree->left->size, tree->left->sizeexp, tree->right->size, tree->right->sizeexp, tree->thread_id);
                } else
                    fprintf(output, "%d, %d, %d, \"int\", %d);\n", tree->left->size, tree->right->size, tree->size, tree->thread_id);
            } else if ((tree->left->arraytype == 1 && tree->right->arraytype == 1) && tree->opid != BOP_dot) {
                if (tree->right->ftype == 1 || tree->left->ftype == 1)
                    fprintf(output, "%d, %d, %d, %d, ", tree->left->size, tree->left->sizeexp, tree->right->size, tree->right->sizeexp);
                else
                    fprintf(output, "%d, %d, ", tree->left->size, tree->right->size);
            }
        }
        break;
    case PREOP:
        if (tree->flag == PUB) {
            fprintf(output, "%s", UOP_symbols[tree->opid]);
            ast_expr_show(tree->left);
        } else {
            arg_str = Str("");
            ast_expr_print(arg_str, tree->left);
            if (tree->opid == UOP_inc)
                fprintf(output, "ss_add_ui(%s, %s, 1)", str_string(arg_str), str_string(arg_str));
            else if (tree->opid == UOP_dec)
                fprintf(output, "ss_sub_ui(%s, %s, 1)", str_string(arg_str), str_string(arg_str));
            str_free(arg_str);
        }
        break;
    case POSTOP:
        if (tree->flag == PUB) {
            ast_expr_show(tree->left);
            fprintf(output, "%s", UOP_symbols[tree->opid]);
        } else {
            arg_str = Str("");
            ast_expr_print(arg_str, tree->left);
            if (tree->opid == UOP_inc)
                fprintf(output, "ss_add_ui(%s, %s, 1)", str_string(arg_str), str_string(arg_str));
            else if (tree->opid == UOP_dec)
                fprintf(output, "ss_sub_ui(%s, %s, 1)", str_string(arg_str), str_string(arg_str));
            str_free(arg_str);
        }
        break;
    case ASS:
        /* for simpile pri to pri assignment */
        if (tree->left->flag == PRI && tree->right->flag == PRI) {
            if (tree->right->index == 0) {
                if (tree->right->type == FUNCCALL) {
                    arg_str = Str("");
                    ast_expr_print(arg_str, tree->right->left);
                    // for compiler defined functions inv and bits
                    if (!strcmp(str_string(arg_str), "inv"))
                        fprintf(output, "%s", "__s->smc_inv");
                    else if (!strcmp(str_string(arg_str), "bits"))
                        fprintf(output, "%s", "__s->smc_bits");
                    else
                        ast_expr_show(tree->right->left);
                    fprintf(output, "(");

                    if (tree->right->right) {
                        ast_expr_show(tree->right->right);
                        fprintf(output, ", ");
                    }

                    ast_expr_show(tree->left);
                    fprintf(output, ")");
                    str_free(arg_str);
                } else {
                    ast_priv_assignment_show(tree, -1);
                }
            }
            // right side is an arithmetic expression
            else
                ast_priv_assignment_show(tree, -1);
            break;
        }
        // left side is private while the right side is public
        else if (tree->left->flag == PRI && tree->right->flag == PUB) {
            if (is_private_struct(tree->left) || is_private_struct_field(tree->left) || is_private_indexed(tree->left)) {
                ast_priv_assignment_show(tree, -1);
                break;
            }
            if (is_priv == 1 && gf == 1) {
                fprintf(global_stream, "__s->smc_set(");
                is_smc_set = 1;
            } else {
                fprintf(output, "__s->smc_set("); // This is where the smc gets printed for private global variables 
            }
            ast_expr_show(tree->right); // this is where the name and value gets printed after smc_set()
            if (is_priv == 1 && gf == 1) {
                fprintf(global_stream, ", ");
            } else {
                fprintf(output, ", ");
            }
            ast_expr_show(tree->left); // this is where the name and value gets printed after smc_set() for op and values I think!

            /* print the bitlength of parameter */
            if (tree->left->ftype == 1) {
                if (is_priv == 1 && gf == 1) {
                    fprintf(global_stream, ", %d, %d, %d, %d, \"float\", %d)", tree->right->size, tree->right->sizeexp, tree->left->size, tree->left->sizeexp, tree->thread_id); // Print it to the global_stream 
                    // printf("%d, %d, %d, %d, \"float\", %d\n", tree->left->size, tree->left->size, tree->left->thread_id); // Print it to the terminal 
                } else {
                    fprintf(output, ", %d, %d, %d, %d, \"float\", %d)", tree->right->size, tree->right->sizeexp, tree->left->size, tree->left->sizeexp, tree->thread_id); // Print to output stream
                }
            } else if (tree->left->ftype == 0) {
                if (is_priv == 1 && gf == 1) {
                    fprintf(global_stream, ", %d, %d, \"int\", %d)", tree->right->size, tree->left->size, tree->thread_id); // Print it to the global_stream 
                    // printf("%d, %d, \"int\", %d\n", tree->left->size, tree->left->size, tree->left->thread_id); // Print it to the terminal 
                } else {
                    fprintf(output, ", %d, %d, \"int\", %d)", tree->right->size, tree->left->size, tree->thread_id); // Print to output stream
                }
            }
            break;
        } else {
            if (tree->right->type == EPMALLOC) {
                ast_expr_show(tree->right);
                if (tree->left->flag == PUB) {
                    ast_expr_show(tree->left);
                    fprintf(output, " %s (struct %s*)_picco_temp_", ASS_symbols[tree->opid], tree->right->u.dtype->spec->name->name);
                } else {
                    char *type = (char *)malloc(sizeof(char) * buffer_size);
                    if (tree->left->u.sym->type == 2)
                        sprintf(type, "struct");
                    else {
                        if (tree->left->ftype == 0)
                            sprintf(type, "int");
                        if (tree->left->ftype == 1)
                            sprintf(type, "float");
                    }
                    fprintf(output, "__s->smc_set_%s_ptr(", type);
                    ast_expr_show(tree->left);
                    fprintf(output, ", (struct %s*)_picco_temp_, \"%s\", %d)", tree->right->u.dtype->spec->name->name, type, tree->right->thread_id);
                    free(type);
                }
            } else {
                // determine if the left operator is private pointer to struct type (consider both pointer variable and field)
                if (is_private_struct(tree->left) || is_private_struct_field(tree->left)) {
                    ast_priv_assignment_show(tree, -1);
                } else {
                    /*for pub assignment */
                    ast_expr_show(tree->left);
                    fprintf(output, " %s ", ASS_symbols[tree->opid]);
                    ast_expr_show(tree->right);
                }
            }
            break;
        }
    case DESIGNATED:
        ast_expr_show(tree->left);
        fprintf(output, " = ");
        ast_expr_show(tree->right);
        break;
    case IDXDES:
        fprintf(output, "[");
        ast_expr_show(tree->left);
        fprintf(output, "]");
        break;
    case DOTDES:
        fprintf(output, ".%s", tree->u.sym->name);
        break;
    case COMMALIST:
    case SPACELIST:
        ast_expr_show(tree->left);
        fprintf(output, "%s", tree->type == COMMALIST ? ", " : " ");
        ast_expr_show(tree->right);
        break;
    case EPMALLOC:
        ast_expr_pmalloc_show(tree);
        break;
    default:
        fprintf(stderr, "[ast_expr_show]: b u g !!\n");
    }
}

void ast_priv_single_expr_show(astexpr e1) {
    if (e1->index > 0) {
        ast_expr_show(e1);
        indent();
    }
    astexpr e0 = Constant(strdup("0"));
    e0->flag = PUB;
    astexpr tree = BinaryOperator(BOP_add, e1, e0);
    tree->ftype = e1->ftype;
    tree->thread_id = e1->thread_id;

    ast_expr_priv_index(tree);
    ast_expr_ptr_dereference(tree, -1);
    ast_expr_refer_struct_field(tree, -1);
}

void ast_priv_cast_helper_show(astexpr tree) {
    if (is_private_indexed(tree) || is_ptr_dereferenced(tree) || is_private_struct_field(tree))
        ast_print_pi_ptr_operator(tree, 1);
    else {
        if (tree->index <= 0)
            ast_expr_show(tree);
        else {
            if (tree->ftype == 1)
                fprintf(output, "_picco_ftmp%d", tree->index);
            else if (tree->ftype == 0)
                fprintf(output, "_picco_tmp%d", tree->index);
        }
        fprintf(output, ", ");
    }
}

// Print specification of tree representing a type or a declaration
void ast_spec_show(astspec tree) {
    switch (tree->type) {
    case SPEC:
        if (tree->subtype == SPEC_int || tree->subtype == SPEC_float)
            fprintf(output, "priv_int");
        else
            fprintf(output, "%s", SPEC_symbols[tree->subtype]);
        break;
    case STCLASSSPEC:
        fprintf(output, "%s", SPEC_symbols[tree->subtype]);
        break;
    case USERTYPE:
        fprintf(output, "%s", tree->name->name);
        break;
    case SUE:
        switch (tree->subtype) {
        case SPEC_enum:
            fprintf(output, "enum");
            if (tree->name)
                fprintf(output, " %s", tree->name->name);
            if (tree->body) {
                fprintf(output, " {\n");
                indlev += 2;
                indent();
                ast_spec_show(tree->body);
                fprintf(output, "\n");
                indlev--;
                indent();
                indlev--;
                fprintf(output, "}");
            }
            break;
        case SPEC_struct:
        case SPEC_union:
            fprintf(output, "%s", tree->subtype == SPEC_struct ? "struct" : "union");
            if (tree->name)
                fprintf(output, " %s", tree->name->name);
            if (tree->u.decl) {
                fprintf(output, " {\n");
                indlev += 2;
                indent();
                tree->u.decl->is_decl_for_sng = tree->is_spec_for_sng;
                ast_decl_show(tree->u.decl);
                fprintf(output, "\n");
                indlev--;
                indent();
                indlev--;
                fprintf(output, "}");
            }
            break;
        default:
            fprintf(stderr, "[ast_spec_show]: SUE b u g !!\n");
        }
        break;
    case ENUMERATOR:
        fprintf(output, "%s", tree->name->name);
        if (tree->u.expr) {
            fprintf(output, " = ");
            ast_expr_show(tree->u.expr);
        }
        break;
    case SPECLIST:
        switch (tree->subtype) {
        case SPEC_Rlist:
            // for "public int"
            if (tree->body->subtype == SPEC_public) {
                fprintf(output, "%s", SPEC_symbols[tree->u.next->subtype]);
                break;
            }
            // for rest of types
            else if (tree->body->subtype != SPEC_private && tree->body->subtype != SPEC_public)
                ast_spec_show(tree->body);
            if (tree->body->type != SPEC || tree->body->subtype != SPEC_star)
                fprintf(output, " ");
            ast_spec_show(tree->u.next);
            break;
        case SPEC_Llist:
            ast_spec_show(tree->u.next);
            fprintf(output, " ");
            ast_spec_show(tree->body);
            break;
        case SPEC_enum:
            ast_spec_show(tree->u.next);
            fprintf(output, ", ");
            ast_spec_show(tree->body);
            break;
        default:
            fprintf(stderr, "[ast_spec_show]: list b u g !!\n");
        }
        break;
    default:
        fprintf(stderr, "[ast_spec_show]: b u g !!\n");
    }
    fflush(stdout);
}

void ast_handle_memory_for_private_variable(astdecl tree, astspec spec, char *struct_name, int flag) {
    switch (tree->type) {
    case DSTRUCTFIELD:
        ast_handle_memory_for_private_variable(tree->decl, tree->spec, struct_name, flag);
        break;
    case DECLARATOR:
        if (tree->decl->type == DARRAY) {
            ast_handle_memory_for_private_variable(tree->decl, spec, struct_name, flag);
            break;
        }
        if (spec->subtype == SPEC_int || spec->subtype == SPEC_Rlist && (spec->body->subtype == SPEC_private && spec->u.next->subtype == SPEC_int)) {
            indent();
            if (tree->spec) {
                if (!flag) {
                    int level = ast_compute_ptr_level(tree);
                    fprintf(output, "%s%s = __s->smc_new_ptr(%d, 0);\n", struct_name, tree->decl->u.id->name, level);
                } else
                    fprintf(output, "__s->smc_free_ptr(&(%s%s));\n", struct_name, tree->decl->u.id->name);
            } else if (!flag)
                fprintf(output, "ss_init(%s%s);\n", struct_name, tree->decl->u.id->name);
            else
                fprintf(output, "ss_clear(%s%s);\n", struct_name, tree->decl->u.id->name);
            break;
        }
        if (spec->subtype == SPEC_float || spec->subtype == SPEC_Rlist && (spec->body->subtype == SPEC_private && spec->u.next->subtype == SPEC_float)) {
            indent();
            if (tree->spec) {
                if (!flag) {
                    int level = ast_compute_ptr_level(tree);
                    fprintf(output, "%s%s = __s->smc_new_ptr(%d, 1);\n", struct_name, tree->decl->u.id->name, level);
                } else
                    fprintf(output, "__s->smc_free_ptr(&(%s%s));\n", struct_name, tree->decl->u.id->name);
            } else {
                if (!flag) {
                    fprintf(output, "%s%s = (priv_int*)malloc(sizeof(priv_int) * (4));\n", struct_name, tree->decl->u.id->name);
                    indent();
                    fprintf(output, "for (int _picco_i = 0; _picco_i < 4; _picco_i++)\n");
                    indlev++;
                    indent();
                    fprintf(output, "ss_init(%s%s[_picco_i]);\n", struct_name, tree->decl->u.id->name);
                    indlev--;
                } else {
                    indent();
                    fprintf(output, "for (int _picco_i = 0; _picco_i < 4; _picco_i++)\n");
                    indlev++;
                    indent();
                    indent();
                    fprintf(output, "ss_clear(%s%s[_picco_i]);\n", struct_name, tree->decl->u.id->name);
                    indlev--;
                    indent();
                    fprintf(output, "free(%s%s);\n", struct_name, tree->decl->u.id->name);
                }
            }
            break;
        }
        if (spec->subtype == SPEC_struct || spec->subtype == SPEC_union) {
            if (!tree->spec) {
                indent();
                if (!flag)
                    fprintf(output, "%s_init(&(%s%s));\n", spec->name->name, struct_name, tree->decl->u.id->name);
                else
                    fprintf(output, "%s_free(&(%s%s));\n", spec->name->name, struct_name, tree->decl->u.id->name);
            } else {
                if (!struct_node_get_flag(sns, spec->name->name)) {
                    indent();
                    if (!flag) {
                        int level = ast_compute_ptr_level(tree);
                        fprintf(output, "%s%s = __s->smc_new_ptr(%d, 2);\n", struct_name, tree->decl->u.id->name, level);
                    } else
                        fprintf(output, "__s->smc_clear_ptr(&(%s%s));\n", struct_name, tree->decl->u.id->name);
                }
            }
            break;
        }
        break;
    case DLIST:
        ast_handle_memory_for_private_variable(tree->u.next, spec, struct_name, flag);
        ast_handle_memory_for_private_variable(tree->decl, spec, struct_name, flag);
        break;
    case DARRAY:
        if (spec->subtype == SPEC_int || spec->subtype == SPEC_Rlist && (spec->body->subtype == SPEC_private && spec->u.next->subtype == SPEC_int)) {
            if (!flag) {
                ast_decl_memory_assign_int(tree, struct_name);
                indlev--;
            } else
                ast_decl_memory_free_int(tree, struct_name);
        } else if (spec->subtype == SPEC_float || spec->subtype == SPEC_Rlist && (spec->body->subtype == SPEC_private && spec->u.next->subtype == SPEC_float)) {
            if (!flag) {
                ast_decl_memory_assign_float(tree, struct_name);
                indlev--;
            } else
                ast_decl_memory_free_float(tree, struct_name);
        }
        break;
    }
}

void ast_return_struct_field_info(char *struct_name, char *field_name, int *ispointer, int *field_type) {
    struct_node node = struct_node_lookup(sns, struct_name);
    struct_field field = struct_field_lookup(node, field_name);
    astdecl f_name = field->name;
    astspec f_type = field->type;
    if (f_name->spec)
        *ispointer = 1;
    else
        *ispointer = 0;
    if (f_type->subtype == SPEC_int || f_type->subtype == SPEC_Rlist && f_type->u.next->subtype == SPEC_int)
        *field_type = 0;
    else if (f_type->subtype == SPEC_float || f_type->subtype == SPEC_Rlist && f_type->u.next->subtype == SPEC_float)
        *field_type = 1;
    else
        *field_type = 2;
}

void ast_print_struct_helper_function(astspec tree) {
    char *init_func = (char *)malloc(sizeof(char) * buffer_size);
    char *free_func = (char *)malloc(sizeof(char) * buffer_size);
    char *type = (char *)malloc(sizeof(char) * buffer_size);
    char *struct_name = (char *)malloc(sizeof(char) * buffer_size);

    if (tree->name) {
        sprintf(init_func, "%s_init", tree->name->name);
        sprintf(free_func, "%s_free", tree->name->name);
        sprintf(struct_name, "%s_node->", tree->name->name);
    }

    if (tree->subtype == SPEC_struct)
        sprintf(type, "struct");
    else
        sprintf(type, "union");

    fprintf(output, "\n");
    fprintf(output, "void %s(%s %s *%s_node)\n", init_func, type, tree->name->name, tree->name->name);
    fprintf(output, "{\n");
    indlev++;
    if (tree->u.decl)
        ast_handle_memory_for_private_variable(tree->u.decl, NULL, struct_name, 0);
    indlev--;
    fprintf(output, "}\n\n");
    fprintf(output, "void %s(%s %s *%s_node)\n", free_func, type, tree->name->name, tree->name->name);
    fprintf(output, "{\n");
    indlev++;
    if (tree->u.decl)
        ast_handle_memory_for_private_variable(tree->u.decl, NULL, struct_name, 1);
    indlev--;
    fprintf(output, "}\n\n");

    // struct field dereference functions
    struct_field_stack sfs = sns->head->fieldlist;
    struct_field field = sfs->head;
    // iterate through the field list of struct
    if (!struct_node_get_flag(sns, tree->name->name)) {
        while (field != NULL) {
            char *field_func_signature = (char *)malloc(sizeof(char) * buffer_size);
            char *field_type_signature = (char *)malloc(sizeof(char) * buffer_size);
            char *struct_ptr = (char *)malloc(sizeof(char) * buffer_size);
            char *var_name = (char *)malloc(sizeof(char) * buffer_size);
            char *array_index = (char *)malloc(sizeof(char) * buffer_size);
            char *field_element = (char *)malloc(sizeof(char) * buffer_size);
            int array_dimension = 0;
            astdecl field_name = field->name;
            astspec field_type = field->type;
            int ispointer = 0;
            int spec_type = 0; // 0-int, 1-float, 2-struct
            int level = 0;
            if (field_name->spec)
                ispointer = 1;
            if (field_type->subtype == SPEC_int || field_type->subtype == SPEC_Rlist && field_type->u.next->subtype == SPEC_int)
                spec_type = 0;
            else if (field_type->subtype == SPEC_float || field_type->subtype == SPEC_Rlist && field_type->u.next->subtype == SPEC_float)
                spec_type = 1;
            else
                spec_type = 2;

            // consider if the struct field is a static array
            // one dimension

            if (field_name->decl->type == DARRAY && field_name->decl->decl->type == DIDENT) {
                array_dimension = 1;
                sprintf(var_name, "%s", field_name->decl->decl->u.id->name);
                sprintf(array_index, "int _picco_index1, ");
                sprintf(field_element, "%s[_picco_index1]", var_name);
            }
            // two-dimension
            else if (field_name->decl->type == DARRAY && field_name->decl->decl->type == DARRAY) {
                array_dimension = 2;
                sprintf(var_name, "%s", field_name->decl->decl->decl->u.id->name);
                sprintf(array_index, "int _picco_index1, int _picco_index2, ");
                sprintf(field_element, "%s[_picco_index1][_picco_index2]", var_name);
            }
            // other cases
            else {
                sprintf(var_name, "%s", field_name->decl->u.id->name);
                sprintf(array_index, "");
                sprintf(field_element, "%s", var_name);
            }

            sprintf(field_func_signature, "%s_%s", tree->name->name, var_name);

            if (!ispointer) {
                if (spec_type == 0)
                    sprintf(field_type_signature, "%s", "priv_int");
                else if (spec_type == 1)
                    sprintf(field_type_signature, "%s", "priv_int*");
            } else {
                sprintf(field_type_signature, "%s", "priv_ptr");
                level = ast_compute_ptr_level(field_name);
            }

            sprintf(struct_ptr, "%s", tree->name->name);
            fprintf(output, "void %s(priv_ptr %s, %s _picco_val, %sint _picco_tag, priv_int _picco_priv_cond, int _picco_thread_id)\n", field_func_signature, struct_ptr, field_type_signature, array_index);
            // print the body of function
            fprintf(output, "{\n");
            indlev++;

            // declare a temporary pointer variable
            char *temp_type = (char *)malloc(sizeof(char) * buffer_size);
            if (spec_type == 0)
                sprintf(temp_type, "int");
            else if (spec_type == 1)
                sprintf(temp_type, "float");
            else if (spec_type == 2)
                sprintf(temp_type, "struct");

            indent();
            fprintf(output, "priv_ptr _picco_tmp_ptr = __s->smc_new_ptr(%d, %d);\n", level + 1, spec_type);
            indent();
            fprintf(output, "listnode _picco_listnode = %s->list->head->next;\n", struct_ptr);
            indent();
            fprintf(output, "while(_picco_listnode != %s->list->tail)\n", struct_ptr);
            indent();
            fprintf(output, "{\n");
            indlev++;

            if (ispointer) {
                indent();
                fprintf(output, "__s->smc_add_ptr(_picco_tmp_ptr, &(((struct %s*)(_picco_listnode->u.struct_var_location))->%s), _picco_listnode->priv_tag, _picco_thread_id);\n", tree->name->name, field_element);
            } else {
                indent();
                fprintf(output, "__s->smc_add_%s_ptr(_picco_tmp_ptr, &(((struct %s*)(_picco_listnode->u.struct_var_location))->%s), _picco_listnode->priv_tag, _picco_thread_id);\n", temp_type, tree->name->name, field_element);
            }

            indent();
            fprintf(output, "_picco_listnode = _picco_listnode->next;\n");
            indlev--;
            indent();
            fprintf(output, "}\n");
            indent();
            fprintf(output, "if(_picco_tag == 0)\n");
            indent();
            fprintf(output, "{\n");
            indlev++;
            if (level > 0) {
                indent();
                fprintf(output, "_picco_val->level = %d;\n", level);
            }
            indent();
            fprintf(output, "__s->smc_dereference_read_ptr(_picco_tmp_ptr, _picco_val, 1, _picco_priv_cond, \"%s\", _picco_thread_id);\n", temp_type);
            indlev--;
            indent();
            fprintf(output, "}\n");
            indent();
            fprintf(output, "else\n");
            indlev++;
            indent();
            fprintf(output, "__s->smc_dereference_write_ptr(_picco_tmp_ptr, _picco_val, 1, _picco_priv_cond, \"%s\", _picco_thread_id);\n", temp_type);
            indlev--;
            indent();
            fprintf(output, "__s->smc_free_ptr(&_picco_tmp_ptr);\n");
            indlev--;
            fprintf(output, "}\n\n");

            field = field->next;

            free(temp_type);
            free(field_func_signature);
            free(field_type_signature);
            free(var_name);
            free(array_index);
            free(field_element);
        }
    }
    free(init_func);
    free(free_func);
    free(type);
    free(struct_name);
}

void ast_decl_sng_stmt_show(aststmt tree) {
    if (tree->is_stmt_for_sng == 2) {
        if (tree->u.declaration.spec->subtype == SPEC_int || tree->u.declaration.spec->subtype == SPEC_float) {
            ast_priv_decl_sng_show(tree->u.declaration.decl, tree->u.declaration.spec);
            return;
        }
        /* for explicit "private" */
        else if (tree->u.declaration.spec->subtype == SPEC_Rlist && tree->u.declaration.spec->body->subtype == SPEC_private) {
            ast_priv_decl_sng_show(tree->u.declaration.decl, tree->u.declaration.spec);
            return;
        } else if ((tree->u.declaration.spec->subtype == SPEC_struct || tree->u.declaration.spec->subtype == SPEC_union)) {
            /* struct variable declaration */
            if (tree->u.declaration.decl) {
                ast_priv_decl_sng_show(tree->u.declaration.decl, tree->u.declaration.spec);
                return;
            }
            /* struct declaration */
            else {
                printf("It is not allowed to declare a struct within an OpenMP construct.\n");
                exit(0);
            }
        }
    }
    tree->u.declaration.spec->is_spec_for_sng = 1;
    ast_spec_show(tree->u.declaration.spec);
    if (tree->u.declaration.decl) {
        fprintf(output, " ");
        ast_decl_show(tree->u.declaration.decl);
    }
    fprintf(output, ";\n");
    return;
}

/**
 * Display of declaration statements in the abstract syntax tree and handle output streams.
 * - Determines if the declaration specifies a global flag.
 * - Handles declarations of int and float types, including private variables.
 * - Processes explicit "private" declarations and updates private variable status.
 * - Handles struct variable and struct declaration cases separately.
 * - Processes explicit "public" and other types of declarations.
 *
 * @param tree The AST representing declaration statements to display.
 * @param current The current branch node in the AST.
 */
void ast_decl_stmt_show(aststmt tree, branchnode current) {
    if (tree->gflag == 1){
        gf = tree->gflag;
    }
    if (tree->u.declaration.spec->subtype == SPEC_int || tree->u.declaration.spec->subtype == SPEC_float)
        ast_priv_decl_show(tree->u.declaration.decl, tree->u.declaration.spec, current, tree->gflag);
    /* for explicit "private" */
    else if (tree->u.declaration.spec->subtype == SPEC_Rlist && tree->u.declaration.spec->body->subtype == SPEC_private) {
        is_priv = 1; // This is where private varibles are known
        ast_priv_decl_show(tree->u.declaration.decl, tree->u.declaration.spec, current, tree->gflag);
    } else if (tree->u.declaration.spec->subtype == SPEC_struct || tree->u.declaration.spec->subtype == SPEC_union) {
        // struct variable declaration
        if (tree->u.declaration.decl) {
            ast_priv_decl_show(tree->u.declaration.decl, tree->u.declaration.spec, current, tree->gflag); // gives error in here 
        }
        // struct declaration
        else {
            struct_node_push(sns, tree->u.declaration.spec);
            struct_node_update(sns, tree->u.declaration.spec->u.decl);
            ast_spec_show(tree->u.declaration.spec);
            fprintf(output, ";\n");
            ast_print_struct_helper_function(tree->u.declaration.spec);
        }
    }
    /* for explicit "public" and others */ 
    else {
        ast_spec_show(tree->u.declaration.spec);
        is_priv = 0;
        if (tree->u.declaration.decl) {
            fprintf(output, " ");
            gf = tree->gflag;
            ast_decl_show(tree->u.declaration.decl);
        }
        fprintf(output, ";\n");
    }
    return;
}

void ast_comma_expr_show(astexpr tree) {
    if (tree->type == COMMALIST) {
        while (1) {
            indent();
            ast_expr_show(tree->right);
            fprintf(output, ";\n");
            if (tree->left->type != COMMALIST)
                break;
            tree = tree->left;
        }

        indent();
        ast_expr_show(tree->left);
        fprintf(output, ";\n");
    } else {
        // indent();
        ast_expr_show(tree);
        fprintf(output, ";\n");
    }
    return;
}

int ast_check_priv_if(astexpr tree) {
    int flag = 0;
    if (tree->type == COMMALIST) {
        while (1) {
            if (tree->right->flag == PRI)
                flag = 1;
            if (tree->left->type != COMMALIST)
                break;
            tree = tree->left;
        }
        if (tree->left->flag == PRI)
            flag = 1;
    } else if (tree->flag == PRI)
        flag = 1;

    return flag;
}

/*
void ast_ptr_decl_show(astdecl tree, astspec spec)
{
    arg_str = Str("");
        ast_expr_print(arg_str, tree->u.expr);
    astdecl d = tree->decl;
    int level = ast_compute_ptr_level(d);
    if(spec->subtype == SPEC_int)
           fprintf(output, "priv_int_ptr %s = __s->smc_new_ptr(%d, 0);\n", tree->decl->decl->u.id->name, level);
        if(spec->subtype == SPEC_float)
           fprintf(output, "priv_float_ptr %s = __s->smc_new_ptr(%d, 1);\n", tree->decl->decl->u.id->name, level);
        indent();

       fprintf(output, "__s->smc_set_ptr(%s, %s, ", tree->decl->decl->u.id->name, str_string(arg_str));
        if(spec->subtype == SPEC_int)
           fprintf(output, "\"int\");\n");
        if(spec->subtype == SPEC_float)
           fprintf(output, "\"float\");\n");
}
*/

int ast_compute_ptr_level(astdecl tree) {
    int level = 1;
    astdecl tmp = tree;
    astspec spec = tree->spec;
    while (spec->type == SPECLIST && spec->body->subtype == SPEC_star) {
        level++;
        spec = spec->u.next;
    }
    return level;
}

void ast_print_priv_field(astdecl tree, astspec spec) {
    switch (tree->type) {
    case DECLARATOR:
        if (tree->decl->type == DARRAY) {
            ast_print_priv_field(tree->decl, spec);
            break;
        }
        if (spec->subtype == SPEC_int) {
            fprintf(output, "priv_int %s", tree->decl->u.id->name);
            break;
        }
        if (spec->subtype == SPEC_float) {
            fprintf(output, "priv_int* %s", tree->decl->u.id->name);
            break;
        }
    case DARRAY:
        if (spec->subtype == SPEC_int) {
            if (tree->decl->type == DIDENT)
                fprintf(output, "priv_int* %s", tree->decl->u.id->name);
            else if (tree->decl->type == DARRAY) {
                if (tree->decl->decl->type == DIDENT)
                    fprintf(output, "priv_int** %s", tree->decl->decl->u.id->name);
            }
        } else if (spec->subtype == SPEC_float) {
            if (tree->decl->type == DIDENT)
                fprintf(output, "priv_int** %s", tree->decl->u.id->name);
            else if (tree->decl->type == DARRAY) {
                if (tree->decl->decl->type == DIDENT)
                    fprintf(output, "priv_int*** %s", tree->decl->decl->u.id->name);
            }
        }
        break;
    }
}

void ast_priv_decl_sng_show(astdecl tree, astspec spec) {
    switch (tree->type) {
    case DECLARATOR:
        if (tree->decl->type == DARRAY) {
            ast_priv_decl_sng_show(tree->decl, spec);
            break;
        }
        if (spec->subtype == SPEC_int || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_int) {
            if (tree->spec)
                fprintf(output, "priv_ptr %s", tree->decl->u.id->name);
            else
                fprintf(output, "priv_int %s", tree->decl->u.id->name);
            break;
        }
        if (spec->subtype == SPEC_float || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_float) {
            if (tree->spec)
                fprintf(output, "priv_ptr %s", tree->decl->u.id->name);
            else
                fprintf(output, "priv_int* %s", tree->decl->u.id->name);
            break;
        }
        /* struct or union declaration -- for now we only support static declaration */
        if (spec->subtype == SPEC_struct || spec->subtype == SPEC_union) {
            /* non-pointer declaration */
            if (!tree->spec) {
                ast_spec_show(spec);
                fprintf(output, " ");
                ast_decl_show(tree);
            }
            /* pointer declaration */
            else {
                /* perform operations below only if the struct contains non-public fields (nested) */
                int contain_pub_field = 1;
                contain_pub_field = struct_node_get_flag(sns, spec->name->name);
                if (!contain_pub_field)
                    fprintf(output, "priv_ptr %s", tree->decl->u.id->name);
                else {
                    ast_spec_show(spec);
                    fprintf(output, " ");
                    ast_decl_show(tree);
                }
            }
            break;
        }
    case DINIT: {
        ast_priv_decl_sng_show(tree->decl, spec);
        fprintf(output, " = ");
        ast_expr_show(tree->u.expr);
        fprintf(output, ";\n");
        break;
    }
    case DARRAY:
        if (spec->subtype == SPEC_int || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_int) {
            if (tree->decl->type == DIDENT)
                fprintf(output, "priv_int* %s", tree->decl->u.id->name);
            else if (tree->decl->type == DARRAY) {
                if (tree->decl->decl->type == DIDENT)
                    fprintf(output, "priv_int** %s", tree->decl->decl->u.id->name);
                else if (tree->decl->decl->type == DARRAY) {
                    printf("We do not support array dimension larger than two.\n");
                    exit(0);
                }
            }
        } else if (spec->subtype == SPEC_float || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_float) {
            if (tree->decl->type == DIDENT)
                fprintf(output, "priv_int** %s", tree->decl->u.id->name);
            else if (tree->decl->type == DARRAY) {
                if (tree->decl->decl->type == DIDENT)
                    fprintf(output, "priv_int*** %s", tree->decl->decl->u.id->name);
                else if (tree->decl->decl->type == DARRAY) {
                    printf("We do not support array dimension larger than two.\n");
                    exit(0);
                }
            }
        }
        break;
    }
}

/**
 * Display/Traverse private declaration statements in the abstract syntax tree and handle output streams.
 * - Processes each declaration in the AST 'tree' based on its type.
 * - Manages the printing of global private declarations into a string.
 * - Handles various types of declarations including int, float, struct, and union.
 * - Manages memory allocation and initialization for private variables.
 * - Supports private arrays with up to two dimensions.
 * - Manages the setup for global flags and private variable status.
 * - Determines whether to print declarations as global or local based on the global flag.
 * - Handles initialization of private variables, including pointers and structs.
 *
 * @param tree The AST representing declaration statements to display.
 * @param spec The AST representing declaration specifications.
 * @param current The current branch node in the AST.
 * @param gflag The global flag indicating whether the declaration is global. (1-means global)
 */
/* Handles private declarations*/
void ast_priv_decl_show(astdecl tree, astspec spec, branchnode current, int gflag) {
    /* printing global private declaration into a string */
    switch (tree->type) {
    case DECLARATOR:
        if(gflag != 1) {
            ltable_push(spec, tree, current->tablelist->head);
        }
        if (tree->decl->type == DARRAY) {
            ast_priv_decl_show(tree->decl, spec, current, gflag);
            break;
        }
        if (spec->subtype == SPEC_int || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_int) {
            // pointer to private int
            if (tree->spec) {
                int level = ast_compute_ptr_level(tree);
                fprintf(output, "priv_ptr %s = __s->smc_new_ptr(%d, 0);\n", tree->decl->u.id->name, level);
            }
            // private int
            else {
                fprintf(output, "priv_int %s;\n", tree->decl->u.id->name); // init in one place, ss_init in other
                indent();
                if (is_priv == 1 && gflag == 1) {
                    fprintf(global_stream, "ss_init(%s);\n", tree->decl->u.id->name);
                } else {
                    fprintf(output, "ss_init(%s);\n", tree->decl->u.id->name);
                }
            }
            break;
        }
        if (spec->subtype == SPEC_float || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_float) {
            // pointer to private float
            if (tree->spec) {
                int level = ast_compute_ptr_level(tree);
                fprintf(output, "priv_ptr %s = __s->smc_new_ptr(%d, 1);\n", tree->decl->u.id->name, level);
            } else { // private float
                fprintf(output, "priv_int* %s; \n", tree->decl->u.id->name); // Float is represented by 4 priv_ints 
                indent();
                // output order for global private 
                // if(gflag != 1){
                if (is_priv == 1 && gflag == 1){
                    fprintf(global_stream, "%s = (priv_int*)malloc(sizeof(priv_int) * (4));\n", tree->decl->u.id->name);
                    indent_global_stream();
                    fprintf(global_stream, "for (int _picco_i = 0; _picco_i < 4; _picco_i++)\n");
                    indlev++;
                    indent_global_stream();
                    indent_global_stream();
                    fprintf(global_stream, "ss_init(%s[_picco_i]);\n", tree->decl->u.id->name);
                    indlev--;
                } else {
                    fprintf(output, "%s = (priv_int*)malloc(sizeof(priv_int) * (4));\n", tree->decl->u.id->name);
                    indent();
                    fprintf(output, "for (int _picco_i = 0; _picco_i < 4; _picco_i++)\n");
                    indlev++;
                    indent();
                    indent();
                    fprintf(output, "ss_init(%s[_picco_i]);\n", tree->decl->u.id->name);
                    indlev--;
                }
            }
            break;
        }
        // struct or union declaration -- for now we only support static declaration
        if (spec->subtype == SPEC_struct || spec->subtype == SPEC_union) {
            // non-pointer declaration
            if (!tree->spec) {
                ast_spec_show(spec);
                fprintf(output, " ");
                ast_decl_show(tree);
                fprintf(output, ";\n");
                indent();
                fprintf(output, "%s_init(&%s);\n", spec->name->name, tree->decl->u.id->name);
            }
            // pointer declaration
            else {
                // perform operations below only if the struct contains non-public fields (nested)
                int contain_pub_field = 1;
                contain_pub_field = struct_node_get_flag(sns, spec->name->name);
                if (!contain_pub_field) {
                    indent();
                    int level = ast_compute_ptr_level(tree);
                    fprintf(output, "priv_ptr %s = __s->smc_new_ptr(%d, 2);\n", tree->decl->u.id->name, level);
                } else {
                    ast_spec_show(spec);
                    fprintf(output, " ");
                    ast_decl_show(tree);
                    fprintf(output, ";\n");
                }
            }
            break;
        }
    case DINIT: {
        ast_priv_decl_show(tree->decl, spec, current, gflag);
        indent();
        astexpr e0 = Identifier(tree->decl->decl->u.id);
        astexpr e1 = Assignment(e0, ASS_eq, tree->u.expr);

        if (spec->subtype == SPEC_int) {
            e0->u.sym->type = 0;
            e0->size = spec->size;
            e0->ftype = e1->ftype = 0;
        }

        if (spec->subtype == SPEC_float) {
            e0->u.sym->type = 1;
            e0->ftype = e1->ftype = 1;
            e0->size = spec->size;
            e0->sizeexp = spec->sizeexp;
        }

        if (spec->subtype == SPEC_struct || spec->subtype == SPEC_union)
            e0->u.sym->type = 2;
        if (spec->subtype == SPEC_Rlist && spec->body->subtype == SPEC_private) {
            if (spec->u.next->subtype == SPEC_int) {
                e0->u.sym->type = 0;
                e0->ftype = e1->ftype = 0;
                e0->size = spec->u.next->size;
            }
            if (spec->u.next->subtype == SPEC_float) {
                e0->u.sym->type = 1;
                e0->ftype = e1->ftype = 1;
                e0->size = spec->u.next->size;
                e0->sizeexp = spec->u.next->sizeexp;
            }
        }
        if (spec->subtype == SPEC_struct || spec->subtype == SPEC_union) {
            int contain_pub_field = 1;
            contain_pub_field = struct_node_get_flag(sns, spec->name->name);
            if (!contain_pub_field)
                e0->flag = e1->flag = PRI;
            else
                e0->flag = e1->flag = PUB;
            e0->isptr = (tree->decl->spec) ? ast_compute_ptr_level(tree->decl) : 0;
            tree->decl->decl->u.id->struct_type = spec;
        } else {
            e0->isptr = (tree->decl->spec) ? ast_compute_ptr_level(tree->decl) : 0;
            e0->flag = e1->flag = PRI;
            // e0->ftype = e1->ftype = tree->u.expr->ftype;
        }
        ast_expr_show(e1);
        if (is_priv == 1 && gflag == 1){
            fprintf(global_stream, ";\n");
        } else {
            fprintf(output, ";\n");
        }
        break;
    }
    // symbol list
    case DLIST:
        ast_priv_decl_show(tree->u.next, spec, current, gflag);
        indent();
        ast_priv_decl_show(tree->decl, spec, current, gflag);
        break;
    // priv array (supports up to two dimensions so far).
    case DARRAY:
        if (spec->subtype == SPEC_int || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_int) {
            if (tree->decl->type == DIDENT)
                fprintf(output, "priv_int* %s; \n", tree->decl->u.id->name);
            // for two-dimensional array
            else if (tree->decl->type == DARRAY) {
                if (tree->decl->decl->type == DIDENT)
                    fprintf(output, "priv_int** %s; \n", tree->decl->decl->u.id->name);
                else if (tree->decl->decl->type == DARRAY) {
                    printf("We do not support array dimension larger than two.\n");
                    exit(0);
                }
            }
            if (symtab_get(stab, tree->decl->u.id, IDNAME) == NULL) {
                ast_decl_memory_assign_int(tree, "");
                indlev--;
            }
        } else if (spec->subtype == SPEC_float || spec->subtype == SPEC_Rlist && spec->u.next->subtype == SPEC_float) {
            if (tree->decl->type == DIDENT)
                fprintf(output, "priv_int** %s; \n", tree->decl->u.id->name);
            else if (tree->decl->type == DARRAY) {
                if (tree->decl->decl->type == DIDENT)
                    fprintf(output, "priv_int*** %s; \n", tree->decl->decl->u.id->name);
                else if (tree->decl->decl->type == DARRAY) {
                    printf("We do not support array dimension larger than two.\n");
                    exit(0);
                }
            }
            if (symtab_get(stab, tree->decl->u.id, IDNAME) == NULL) {
                ast_decl_memory_assign_float(tree, "");
                indlev--;
            }
        }
        break;
    }
}

void ast_decl_memory_assign_int(astdecl tree, char *prefix) {
    indent();
    if (tree->decl->type == DIDENT) {
        fprintf(output, "%s%s = (priv_int*)malloc(sizeof(priv_int) * (", prefix, tree->decl->u.id->name);
        ast_expr_show(tree->u.expr);
    } else if (tree->decl->type == DARRAY) {
        fprintf(output, "%s%s = (priv_int**)malloc(sizeof(priv_int*) * (", prefix, tree->decl->decl->u.id->name);
        ast_expr_show(tree->decl->u.expr);
    }
    fprintf(output, "));\n");
    // init
    indent();
    fprintf(output, "for (int _picco_i = 0; _picco_i < ");
    if (tree->decl->type == DIDENT)
        ast_expr_show(tree->u.expr);
    else if (tree->decl->type == DARRAY)
        ast_expr_show(tree->decl->u.expr);
    fprintf(output, "; _picco_i++)\n");
    indlev++;
    indent();
    if (tree->decl->type == DIDENT) {
        indent();
        fprintf(output, "ss_init(%s%s[_picco_i]);\n", prefix, tree->decl->u.id->name);
    } else if (tree->decl->type == DARRAY) {
        arg_str = Str("");
        ast_expr_print(arg_str, tree->u.expr);
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "%s%s[_picco_i] = (priv_int*)malloc(sizeof(priv_int) * (%s));\n", prefix, tree->decl->decl->u.id->name, str_string(arg_str));
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < %s; _picco_j++)\n", str_string(arg_str));
        indent();
        indent();
        fprintf(output, "ss_init(%s%s[_picco_i][_picco_j]);\n", prefix, tree->decl->decl->u.id->name);
        indlev--;
        indent();
        fprintf(output, "}\n");
        str_free(arg_str);
    }
}

void ast_decl_memory_free_int(astdecl tree, char *prefix) {
    // free
    indent();
    fprintf(output, "for (int _picco_i = 0; _picco_i < ");
    if (tree->decl->type == DIDENT)
        ast_expr_show(tree->u.expr);
    else if (tree->decl->type == DARRAY)
        ast_expr_show(tree->decl->u.expr);
    fprintf(output, "; _picco_i++)\n");
    indlev++;
    indent();
    if (tree->decl->type == DIDENT) {
        indent();
        fprintf(output, "ss_clear(%s%s[_picco_i]);\n", prefix, tree->decl->u.id->name);
        indlev--;
        indent();
        fprintf(output, "free(%s%s);\n", prefix, tree->decl->u.id->name);
    } else if (tree->decl->type == DARRAY) {
        arg_str = Str("");
        ast_expr_print(arg_str, tree->u.expr);
        fprintf(output, "{\n");
        indlev++;
        indent();
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < %s; _picco_j++)\n", str_string(arg_str));
        indent();
        indlev++;
        indent();
        fprintf(output, "ss_clear(%s%s[_picco_i][_picco_j]);\n", prefix, tree->decl->decl->u.id->name);
        indlev--;
        indent();
        indent();
        fprintf(output, "free(%s%s[_picco_i]);\n", prefix, tree->decl->decl->u.id->name, str_string(arg_str));
        indent();
        indlev--;
        fprintf(output, "}\n");
        indlev--;
        indent();
        fprintf(output, "free(%s%s);\n", prefix, tree->decl->decl->u.id->name);
        str_free(arg_str);
    }
}

void ast_decl_memory_assign_float(astdecl tree, char *prefix) {
    indent();
    if (tree->decl->type == DIDENT) {
        fprintf(output, "%s%s = (priv_int**)malloc(sizeof(priv_int*) * (", prefix, tree->decl->u.id->name);
        ast_expr_show(tree->u.expr);
    } else if (tree->decl->type == DARRAY) {
        fprintf(output, "%s%s = (priv_int***)malloc(sizeof(priv_int**) * (", prefix, tree->decl->decl->u.id->name);
        ast_expr_show(tree->decl->u.expr);
    }
    fprintf(output, "));\n");
    // init
    indent();
    fprintf(output, "for (int _picco_i = 0; _picco_i < ");
    if (tree->decl->type == DIDENT)
        ast_expr_show(tree->u.expr);
    else if (tree->decl->type == DARRAY)
        ast_expr_show(tree->decl->u.expr);
    fprintf(output, "; _picco_i++)\n");
    indlev++;
    indent();
    if (tree->decl->type == DIDENT) {
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "%s%s[_picco_i] = (priv_int*)malloc(sizeof(priv_int) * (4));\n", prefix, tree->decl->u.id->name);
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
        indent();
        indent();
        fprintf(output, "ss_init(%s%s[_picco_i][_picco_j]);\n", prefix, tree->decl->u.id->name);
        indlev--;
        indent();
        fprintf(output, "}\n");
    } else if (tree->decl->type == DARRAY) {
        arg_str = Str("");
        ast_expr_print(arg_str, tree->u.expr);
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "%s%s[_picco_i] = (priv_int**)malloc(sizeof(priv_int*) * (%s));\n", prefix, tree->decl->decl->u.id->name, str_string(arg_str));
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < %s; _picco_j++)\n", str_string(arg_str));
        indent();
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "%s%s[_picco_i][_picco_j] = (priv_int*)malloc(sizeof(priv_int) * (4));\n", prefix, tree->decl->decl->u.id->name);
        indent();
        fprintf(output, "for (int _picco_k = 0; _picco_k < 4; _picco_k++)\n");
        indent();
        indent();
        fprintf(output, "ss_init(%s%s[_picco_i][_picco_j][_picco_k]);\n", prefix, tree->decl->decl->u.id->name);
        indlev--;
        indent();
        fprintf(output, "}\n");
        indlev--;
        indent();
        fprintf(output, "}\n");
        str_free(arg_str);
    }
}

void ast_decl_memory_free_float(astdecl tree, char *prefix) {
    // free
    indent();
    fprintf(output, "for (int _picco_i = 0; _picco_i < ");
    if (tree->decl->type == DIDENT)
        ast_expr_show(tree->u.expr);
    else if (tree->decl->type == DARRAY)
        ast_expr_show(tree->decl->u.expr);
    fprintf(output, "; _picco_i++)\n");
    indlev++;
    indent();
    if (tree->decl->type == DIDENT) {
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
        indent();
        indent();
        fprintf(output, "ss_clear(%s%s[_picco_i][_picco_j]);\n", prefix, tree->decl->u.id->name);
        indlev--;
        indent();
        fprintf(output, "free(%s%s[_picco_i]);\n", prefix, tree->decl->u.id->name);
        indent();
        fprintf(output, "}\n");
        indent();
        fprintf(output, "free(%s%s);\n", prefix, tree->decl->u.id->name);
    } else if (tree->decl->type == DARRAY) {
        arg_str = Str("");
        ast_expr_print(arg_str, tree->u.expr);
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "for (int _picco_j = 0; _picco_j < %s; _picco_j++)\n", str_string(arg_str));
        indent();
        fprintf(output, "{\n");
        indlev++;
        indent();
        fprintf(output, "for (int _picco_k = 0; _picco_k < 4; _picco_k++)\n");
        indent();
        indent();
        fprintf(output, "ss_clear(%s%s[_picco_i][_picco_j][_picco_k]);\n", prefix, tree->decl->decl->u.id->name);
        indent();
        fprintf(output, "free(%s%s[_picco_i][_picco_j]);\n", prefix, tree->decl->decl->u.id->name);
        indlev--;
        indent();
        fprintf(output, "}\n");
        indlev--;
        indent();
        fprintf(output, "free(%s%s[_picco_i]);\n", prefix, tree->decl->decl->u.id->name);
        indent();
        fprintf(output, "}\n");
        indent();
        fprintf(output, "free(%s%s);\n", prefix, tree->decl->decl->u.id->name);
        str_free(arg_str);
    }
}

void ast_list_element_show(astdecl tree) {
    if (tree->decl && tree->decl->spec && (tree->spec->subtype != SPEC_struct && tree->spec->subtype != SPEC_union) && (tree->spec->subtype == SPEC_int || tree->spec->subtype == SPEC_float || (tree->spec->subtype == SPEC_Rlist && tree->spec->body->subtype == SPEC_private))) {
        fprintf(output, "priv_ptr ");
        ast_decl_show(tree->decl->decl);
    } else if (tree->decl && tree->decl->spec && tree->spec->subtype == SPEC_struct && !struct_node_get_flag(sns, tree->spec->name->name)) {
        fprintf(output, "priv_ptr ");
        ast_decl_show(tree->decl->decl);
    } else {
        // when the member is a type of private variables and arrays
        if ((tree->spec->subtype == SPEC_int || tree->spec->subtype == SPEC_float) && tree->decl)
            ast_print_priv_field(tree->decl, tree->spec);
        else if (tree->spec->subtype == SPEC_Rlist && tree->spec->body->subtype == SPEC_private && tree->decl)
            ast_print_priv_field(tree->decl, tree->spec->u.next);
        // when the member is a public type
        else {
            ast_spec_show(tree->spec);
            fprintf(output, " ");
            if (tree->decl)
                ast_decl_show(tree->decl);
        }
    }
}

void ast_decl_show(astdecl tree) {
    switch (tree->type) {
    case DIDENT:
        if (!strcmp(tree->u.id->name, "main")) {
            fprintf(output, "%s", " __original_main");
        } else {
            fprintf(output, "%s", tree->u.id->name); // the public variable name gets printed here
        }
        break;
    case DPAREN:
        fprintf(output, "(");
        ast_decl_show(tree->decl);
        fprintf(output, ")");
        break;
    case DARRAY:
        if (tree->decl) /* Maybe abstract declarator */
            ast_decl_show(tree->decl);
        fprintf(output, "[");
        if (tree->spec)
            ast_spec_show(tree->spec);
        if (tree->u.expr) {
            ast_expr_show(tree->u.expr);
        }
        fprintf(output, "]");
        break;
    case DFUNC:
        if (tree->decl)
            ast_decl_show(tree->decl);
        fprintf(output, "(");
        if (tree->u.params)
            ast_decl_show(tree->u.params);
        fprintf(output, ")");
        break;
    case DINIT:
        ast_decl_show(tree->decl);
        if (tree->u.expr != NULL) {
            fprintf(output, " = ");
            ast_expr_show(tree->u.expr);
        }
        break;
    case DECLARATOR:
        if (tree->spec) /* pointer */
        {
            ast_spec_show(tree->spec);
        }
        ast_decl_show(tree->decl);
        break;
    case ABSDECLARATOR:
        if (tree->spec) /* pointer */
            ast_spec_show(tree->spec);
        if (tree->decl) {
            if (tree->spec)
                fprintf(output, " ");
            ast_decl_show(tree->decl);
        }
        break;
    case DPARAM:
        ast_list_element_show(tree);
        break;
    case DELLIPSIS:
        fprintf(output, "...");
        break;
    case DBIT:
        if (tree->decl)
            ast_decl_show(tree->decl);
        fprintf(output, " : ");
        ast_expr_show(tree->u.expr);
        break;
    case DSTRUCTFIELD:
        ast_list_element_show(tree);
        if (!tree->is_decl_for_sng) // if the field does not belong to Openmp struct
            struct_field_push(sns->head->fieldlist, tree->spec, tree->decl);
        fprintf(output, ";");
        break;
    case DCASTTYPE:
        ast_spec_show(tree->spec);
        if (tree->decl) {
            fprintf(output, " ");
            ast_decl_show(tree->decl);
        }
        break;
    case DLIST:
        switch (tree->subtype) {
        case DECL_decllist:
        case DECL_idlist:
        case DECL_paramlist:
            if (tree->u.next == NULL || tree->decl == NULL) {
                fprintf(stderr, "[ast_decl_show]: list next/body NULL !!\n");
                break;
            }
            ast_decl_show(tree->u.next);
            fprintf(output, ", ");
            ast_decl_show(tree->decl);
            break;
        case DECL_fieldlist:
            tree->u.next->is_decl_for_sng = tree->decl->is_decl_for_sng = tree->is_decl_for_sng;
            ast_decl_show(tree->u.next);
            fprintf(output, "\n");
            indent();
            ast_decl_show(tree->decl);
            break;
        default:
            fprintf(stderr, "[ast_decl_show]: list b u g !!\n");
        }
        break;
    default:
        fprintf(stderr, "[ast_decl_show]: b u g !!\n");
    }
    fflush(stdout);
}

void ast_tmp_decl_show(char *prefix, int sidx, int eidx) {
    if (eidx >= sidx) {
        indent();
        fprintf(output, "priv_int ");
        for (ind = sidx; ind < eidx; ind++)
            fprintf(output, "%stmp%d, ", prefix, ind);
        fprintf(output, "%stmp%d;\n", prefix, eidx);
        for (ind = sidx; ind <= eidx; ind++) {
            indent();
            fprintf(output, "ss_init(%stmp%d);\n", prefix, ind);
        }
        fprintf(output, "\n");
    }
}

void ast_float_tmp_decl_show(char *prefix, int sidx, int eidx) {
    if (eidx >= sidx) {
        for (ind = sidx; ind <= eidx; ind++) {
            indent();
            fprintf(output, "priv_int* %stmp%d = (priv_int*)malloc(sizeof(priv_int) * 4);\n", prefix, ind);
            indent();
            fprintf(output, "for(int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
            indlev++;
            indent();
            fprintf(output, "ss_init(%stmp%d[_picco_j]);\n", prefix, ind);
            indlev--;
        }
    }
}

void ast_float_tmp_clear_show(char *prefix, int sidx, int eidx) {
    if (eidx > 0) {
        for (ind = sidx; ind <= eidx; ind++) {
            indent();
            fprintf(output, "for(int _picco_j = 0; _picco_j < 4; _picco_j++)\n");
            indlev++;
            indent();
            fprintf(output, "ss_clear(%s%d[_picco_j]);\n", prefix, ind);
            indlev--;
            indent();
            fprintf(output, "free(%s%d);\n", prefix, ind);
        }
    }
}

void ast_tmp_clear_show(char *prefix, int sidx, int eidx) {
    if (eidx > 0) {
        for (ind = sidx; ind <= eidx; ind++) {
            indent();
            fprintf(output, "ss_clear(%s%d);\n", prefix, ind);
        }
    }
}

void ast_smc_show(char *prefix, astexpr tree) {
    fprintf(output, "%s(", prefix);
}

void ast_priv_assignment_show(astexpr tree, int private_if_index) {
    str rightop = Str("");
    str leftop = Str("");
    char *right = (char *)malloc(sizeof(char) * buffer_size);
    char *left = (char *)malloc(sizeof(char) * buffer_size);
    char *tmp = (char *)malloc(sizeof(char) * buffer_size);
    char *type = (char *)malloc(sizeof(char) * buffer_size);
    char *priv_cond = (char *)malloc(sizeof(char) * buffer_size);
    int dereferences = 0, pointers = 0;
    if (tree->right->index == 0) {
        if (is_private_indexed(tree->right) || is_ptr_dereferenced(tree->right) || is_private_struct_field(tree->right)) {
            astexpr e0 = Constant(strdup("0"));
            astexpr e1 = BinaryOperator(BOP_add, tree->right, e0);
            e1->ftype = tree->right->ftype;
            e1->thread_id = tree->right->thread_id;
            // compute the left operator
            if (is_private_indexed(tree->right))
                ast_expr_priv_index(e1);
            if (is_ptr_dereferenced(tree->right)) {
                ast_expr_ptr_dereference(e1, private_if_index);
                ast_compute_dereferences_and_pointers(e1->left, &dereferences, &pointers);
            }
            if (is_private_struct_field(tree->right)) {
                ast_expr_refer_struct_field(e1, private_if_index);
            }
            // print the left operator
            if (is_private_indexed(tree->right) || is_ptr_dereferenced(tree->right)) {
                if (tree->right->ftype == 0) {
                    if (pointers == dereferences)
                        sprintf(right, "_picco_priv_tmp1");
                    else
                        sprintf(right, "_picco_tmp_int_ptr1");
                } else if (tree->right->ftype == 1) {
                    if (pointers == dereferences)
                        sprintf(right, "_picco_priv_ftmp1");
                    else
                        sprintf(right, "_picco_tmp_float_ptr1");
                }
            }
            // for private struct field
            else if (is_private_struct_field(tree->right)) {
                int ispointer, field_type;
                char *struct_name = (char *)malloc(sizeof(char) * buffer_size);
                astexpr tree1 = NULL;
                if (tree->right->type == ARRAYIDX && tree->right->left->type == PTRFIELD)
                    tree1 = tree->right->left;
                else if (tree->right->type == ARRAYIDX && tree->right->left->type == ARRAYIDX)
                    tree1 = tree->right->left->left;
                else
                    tree1 = tree->right;
                if (tree1->left->type == PTRFIELD || tree1->left->type == DOTFIELD) {
                    astspec type = ast_get_struct_name(tree1->left);
                    sprintf(struct_name, "%s", type->name->name);
                } else {
                    sprintf(struct_name, "%s", tree1->left->u.sym->struct_type->name->name);
                }
                ast_return_struct_field_info(struct_name, tree1->u.sym->name, &ispointer, &field_type);
                ast_declare_temp_for_struct_field(&right, ispointer, field_type, 1);
                free(struct_name);
            }
            astexpr name = String(right);
            ast_expr_print(rightop, name);
        } else {
            ast_expr_print(rightop, tree->right);
        }
    }
    // if right side is an arithmetic operation
    else if (tree->right->index > 0) {
        if (tree->right->ftype == 1)
            sprintf(tmp, "_picco_ftmp%d", tree->right->index);
        else if (tree->right->ftype == 0)
            sprintf(tmp, "_picco_tmp%d", tree->right->index);
        astexpr name = String(tmp);
        ast_expr_print(rightop, name);
        ast_expr_show(tree->right);
    }
    // only corresponds to the case of struct_pointer = 0; - initialization
    else if (tree->right->index < 0) {
        ast_expr_print(rightop, tree->right);
    }
    // considering the left operator
    if (is_private_indexed(tree->left)) {
        char *priv_cond = (char *)malloc(sizeof(char) * buffer_size);
        if (private_if_index != -1)
            sprintf(priv_cond, "_picco_condtmp%d", private_if_index);
        else
            sprintf(priv_cond, "NULL");
        sprintf(left, "_picco_priv_ind3");
        ast_print_private_index(tree->left, left);
        ast_print_private_indexed_stmt_write(tree->left, left, str_string(rightop), -1, tree->opid, priv_cond, NULL, NULL);
        free(priv_cond);
    }
    /***********************************************/
    else if (is_private_struct_field(tree->left)) {
        int recursive_level = 0, tag = 0, id3 = 3;
        int is_address = 0;
        char *var = (char *)malloc(sizeof(char) * buffer_size);
        if (tree->right->type == UOP && tree->right->opid == UOP_addr)
            is_address = 1;
        ast_expr_refer_single_struct_field(tree->left, &id3, is_address, 1, private_if_index, str_string(rightop), &recursive_level, &var, &tag, tree->left->thread_id);
        free(var);
    }
    /**********************************************/
    else if (is_ptr_dereferenced(tree->left)) {
        strip_off_bracket_and_dereferences(leftop, tree->left);
        ast_compute_dereferences_and_pointers(tree->left, &dereferences, &pointers);
        if (private_if_index != -1)
            sprintf(priv_cond, "_picco_condtmp%d", private_if_index);
        else
            sprintf(priv_cond, "NULL");
        indent();
        char *right_type = (char *)malloc(sizeof(char) * buffer_size);
        ast_return_type_of_rop(&right_type, tree->right);
        if (tree->left->ftype == 0)
            fprintf(output, "__s->smc_dereference_write%s_ptr(%s, %s, %d, %s, \"int\", %d)", right_type, str_string(leftop), str_string(rightop), dereferences, priv_cond, tree->left->thread_id);
        else if (tree->left->ftype == 1)
            fprintf(output, "__s->smc_dereference_write%s_ptr(%s, %s, %d, %s, \"float\", %d)", right_type, str_string(leftop), str_string(rightop), dereferences, priv_cond, tree->left->thread_id);
        free(right_type);
    } else if (!is_ptr_dereferenced(tree->left) && tree->left->isptr > 0) {
        ast_ptr_assignment_show(tree->left, tree->right, rightop);
    } else {
        ast_expr_print(leftop, tree->left);
        // perform vector computation
        if (tree->left->arraytype == 1) {
            fprintf(output, "%s, ", str_string(leftop));
            if (tree->right->ftype == 1) {
                if (tree->right->opid != BOP_lt && tree->right->opid != BOP_gt && tree->right->opid != BOP_leq && tree->right->opid != BOP_geq && tree->right->opid != BOP_eqeq && tree->right->opid != BOP_neq)
                    fprintf(output, "%d, %d, ", tree->left->size, tree->left->sizeexp);
                else
                    fprintf(output, "1, ");
            } else
                fprintf(output, "%d, ", tree->left->size);
            ast_expr_show(tree->left->arraysize);

            if (tree->right->left->ftype == 1 || tree->right->right->ftype == 1)
                type = "\"float\"";
            else
                type = "\"int\"";

            fprintf(output, ", %s, %d)", type, tree->right->thread_id);
        } else if (tree->right->opid == BOP_dot) {
            fprintf(output, "%s, %d)", str_string(leftop), tree->right->thread_id);
        } else {
            if (tree->right->ftype == 1)
                type = "\"float\"";
            else if (tree->right->ftype == 0)
                type = "\"int\"";
            indent();
            ast_assignment_prefix_show(tree);
            if (tree->opid != ASS_eq) {
                if (tree->right->ftype == 1)
                    fprintf(output, "(%s, %s, %s, %d, %d, %d, %d, %d, %d, \"float\", %d)", str_string(leftop), str_string(rightop), str_string(leftop), tree->left->size, tree->left->sizeexp, tree->right->size, tree->right->sizeexp, tree->left->size, tree->left->sizeexp, tree->right->thread_id);
                else
                    fprintf(output, "(%s, %s, %s, %d, %d, %d, \"int\", %d)", str_string(leftop), str_string(rightop), str_string(leftop), tree->left->size, tree->right->size, tree->left->size, tree->right->thread_id);
            } else {
                if (tree->right->ftype == 1)
                    fprintf(output, "(%s, %s, %d, %d, %d, %d, \"float\", %d)", str_string(rightop), str_string(leftop), tree->right->size, tree->right->sizeexp, tree->left->size, tree->left->sizeexp, tree->right->thread_id);
                else
                    fprintf(output, "(%s, %s, %d, %d, \"int\", %d)", str_string(rightop), str_string(leftop), tree->right->size, tree->left->size, tree->right->thread_id);
            }
        }
    }
}

void ast_temporary_variable_declaration() {
    fprintf(output, "\n");
    if (tmp_index >= 1)
        ast_tmp_decl_show("_picco_", 1, tmp_index);
    if (tmp_float_index >= 1)
        ast_float_tmp_decl_show("_picco_f", 1, tmp_float_index);
    indent();
    fprintf(output, "void* _picco_temp_;\n ");
    if (is_priv_int_index_appear || is_priv_float_index_appear) {
        indent();
        fprintf(output, "priv_int _picco_priv_ind1, _picco_priv_ind2, _picco_priv_ind3;\n");
        indent();
        fprintf(output, "ss_init(_picco_priv_ind1);\n");
        indent();
        fprintf(output, "ss_init(_picco_priv_ind2);\n");
        indent();
        fprintf(output, "ss_init(_picco_priv_ind3);\n");
    }

    if (is_priv_int_index_appear || is_priv_int_ptr_appear) {
        indent();
        fprintf(output, "priv_int _picco_priv_tmp1, _picco_priv_tmp2;\n");
        indent();
        fprintf(output, "ss_init(_picco_priv_tmp1);\n");
        indent();
        fprintf(output, "ss_init(_picco_priv_tmp2);\n");
    }
    if (is_priv_float_index_appear || is_priv_float_ptr_appear) {
        indent();
        fprintf(output, "priv_int *_picco_priv_ftmp1 = (priv_int*)malloc(sizeof(priv_int) * 4);\n");
        indent();
        fprintf(output, "priv_int *_picco_priv_ftmp2 = (priv_int*)malloc(sizeof(priv_int) * 4);\n");
        indent();
        fprintf(output, "for(int i = 0; i < 4; i++){\n");
        indlev++;
        indent();
        fprintf(output, "ss_init(_picco_priv_ftmp1[i]);\n");
        indent();
        fprintf(output, "ss_init(_picco_priv_ftmp2[i]);\n");
        indlev--;
        fprintf(output, "}\n");
    }
    if (is_priv_int_ptr_appear) {
        indent();
        fprintf(output, "priv_ptr _picco_tmp_int_ptr1 = __s->smc_new_ptr(1, 0);\n");
        indent();
        fprintf(output, "priv_ptr _picco_tmp_int_ptr2 = __s->smc_new_ptr(1, 0);\n");
    }
    if (is_priv_float_ptr_appear) {
        indent();
        fprintf(output, "priv_ptr _picco_tmp_float_ptr1 = __s->smc_new_ptr(1, 1);\n");
        indent();
        fprintf(output, "priv_ptr _picco_tmp_float_ptr2 = __s->smc_new_ptr(1, 1);\n");
    }

    if (is_priv_int_struct_field_appear) {
        indent();
        fprintf(output, "priv_int _picco_str_field_tmp_int1, _picco_str_field_tmp_int2, _picco_str_field_tmp_int3;\n");
        indent();
        fprintf(output, "ss_init(_picco_str_field_tmp_int1);\n");
        indent();
        fprintf(output, "ss_init(_picco_str_field_tmp_int2);\n");
        indent();
        fprintf(output, "ss_init(_picco_str_field_tmp_int3);\n");
    }

    if (is_priv_int_ptr_struct_field_appear) {
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_int_ptr1 = __s->smc_new_ptr(1, 0);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_int_ptr2 = __s->smc_new_ptr(1, 0);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_int_ptr3 = __s->smc_new_ptr(1, 0);\n");
        indent();
        fprintf(output, "priv_ptr _picco_tmp_int_ptr = __s->smc_new_ptr(0, 0);\n");
    }

    if (is_priv_float_struct_field_appear) {

        indent();
        fprintf(output, "priv_int *_picco_str_field_tmp_float1 = (priv_int*)malloc(sizeof(priv_int) * 4);\n");
        indent();
        fprintf(output, "priv_int *_picco_str_field_tmp_float2 = (priv_int*)malloc(sizeof(priv_int) * 4);\n");
        indent();
        fprintf(output, "priv_int *_picco_str_field_tmp_float3 = (priv_int*)malloc(sizeof(priv_int) * 4);\n");
        indent();
        fprintf(output, "for(int i = 0; i < 4; i++){\n");
        indlev++;
        indent();
        fprintf(output, "ss_init(_picco_str_field_tmp_float1[i]);\n");
        indent();
        fprintf(output, "ss_init(_picco_str_field_tmp_float2[i]);\n");
        indent();
        fprintf(output, "ss_init(_picco_str_field_tmp_float3[i]);\n");
        indlev--;
        fprintf(output, "}\n");
    }

    if (is_priv_float_ptr_struct_field_appear) {
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_float_ptr1 = __s->smc_new_ptr(1, 1);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_float_ptr2 = __s->smc_new_ptr(1, 1);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_float_ptr3 = __s->smc_new_ptr(1, 1);\n");
        indent();
        fprintf(output, "priv_ptr _picco_tmp_float_ptr = __s->smc_new_ptr(1, 1);\n");
    }

    if (is_priv_struct_ptr_struct_field_appear) {
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_struct_ptr1 = __s->smc_new_ptr(1, 2);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_struct_ptr2 = __s->smc_new_ptr(1, 2);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_struct_ptr3 = __s->smc_new_ptr(1, 2);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_struct_ptr4 = __s->smc_new_ptr(1, 2);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_struct_ptr5 = __s->smc_new_ptr(1, 2);\n");
        indent();
        fprintf(output, "priv_ptr _picco_str_field_tmp_struct_ptr6 = __s->smc_new_ptr(1, 2);\n");
        indent();
        fprintf(output, "priv_ptr _picco_tmp_struct_ptr = __s->smc_new_ptr(1, 2);\n");
    }
}

void ast_assignment_prefix_show(astexpr tree) {
    switch (tree->opid) {
    case ASS_eq:
        fprintf(output, "__s->smc_set");
        break;
    case ASS_add:
        fprintf(output, "__s->smc_add");
        break;
    case ASS_sub:
        fprintf(output, "__s->smc_sub");
        break;
    case ASS_mul:
        fprintf(output, "__s->smc_mult");
        break;
    case ASS_div:
        fprintf(output, "__s->smc_div");
        break;
    case ASS_mod:
        fprintf(output, "__s->smc_mod");
        break;
    default:
        break;
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                               *
 *     OpenMP NODES                                              *
 *                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void ast_ompclause_show(ompclause t) {
    if (t->type == OCLIST) {
        if (t->u.list.next != NULL) {
            ast_ompclause_show(t->u.list.next);
            fprintf(output, " ");
        }
        assert((t = t->u.list.elem) != NULL);
    }

    fprintf(output, "%s", clausenames[t->type]);
    switch (t->type) {
    case OCIF:
    case OCFINAL:
    case OCNUMTHREADS:
        fprintf(output, "( ");
        ast_expr_show(t->u.expr);
        fprintf(output, " )");
        break;
    case OCSCHEDULE:
        fprintf(output, "( %s%s", clausesubs[t->subtype], t->u.expr ? ", " : " ");
        if (t->u.expr)
            ast_expr_show(t->u.expr);
        fprintf(output, " )");
        break;
    case OCDEFAULT:
        fprintf(output, "( %s )", clausesubs[t->subtype]);
        break;
    case OCREDUCTION:
        fprintf(output, "( %s : ", clausesubs[t->subtype]);
        ast_decl_show(t->u.varlist);
        fprintf(output, ")");
        break;
    case OCCOPYIN:
    case OCPRIVATE:
    case OCCOPYPRIVATE:
    case OCFIRSTPRIVATE:
    case OCLASTPRIVATE:
    case OCSHARED:
        fprintf(output, "(");
        ast_decl_show(t->u.varlist);
        fprintf(output, ")");
        break;
    case OCNOWAIT:
    case OCORDERED:
    case OCUNTIED:
    case OCMERGEABLE:
        break;
    case OCCOLLAPSE:
        fprintf(output, "(%d)", t->subtype);
        break;
    }
}

void ast_ompdir_show(ompdir t) {
    fprintf(output, "#pragma omp %s ", ompdirnames[t->type]);
    switch (t->type) {
    case DCCRITICAL:
        if (t->u.region)
            fprintf(output, "(%s)", t->u.region->name);
        break;
    case DCFLUSH:
        if (t->u.varlist) {
            fprintf(output, "(");
            ast_decl_show(t->u.varlist);
            fprintf(output, ")");
        }
        break;
    case DCTHREADPRIVATE:
        if (t->u.varlist) {
            fprintf(output, "(");
            ast_decl_show(t->u.varlist);
            fprintf(output, ")");
        }
        break;
    default:
        if (t->clauses)
            ast_ompclause_show(t->clauses);
        break;
    }
    fprintf(output, "\n");
}

void ast_ompcon_show(ompcon t, branchnode current) {
    ast_ompdir_show(t->directive);
    if (t->body) /* barrier & flush don't have a body. */
    {
        indent();
        ast_stmt_show(t->body, current);
    }
}

/**
 * Read the content of a file into a dynamically allocated string.
 *
 * This function opens the file "PICCO_Global.txt" for reading. It determines the size of the file,
 * allocates memory to hold its content, reads the file into the allocated memory, and returns a pointer
 * to the allocated string containing the file content. If any errors occur during file operations
 * or memory allocation, the function prints an error message and exits the program.
 *
 * @return A pointer to a dynamically allocated string containing the content of the file.
 */
char *readFileToString() {
    FILE *file = fopen("PICCO_Global.txt", "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    // Seek to end of file to determine file size
    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    fseek(file, 0, SEEK_SET);

    // Allocate memory for the string
    char *fileContent = (char *)malloc(fileSize + 1);
    if (fileContent == NULL) {
        perror("Error allocating memory");
        fclose(file);
        exit(1);
    }

    // Read the file into the string
    size_t bytesRead = fread(fileContent, 1, fileSize, file);
    if (bytesRead != fileSize) {
        perror("Error reading file");
        fclose(file);
        free(fileContent);
        exit(1);
    }
    
    fileContent[bytesRead] = '\0';
    fclose(file);
    
    // Deleting the file
    char remove_file[] = "PICCO_Global.txt";
    if (remove(remove_file) != 0) {
        printf("Error deleting the file %s.\n", remove_file);
    }
    
    return fileContent;
}


/**
 * Display an abstract syntax tree and handle output streams.
 *
 * This function opens the output stream and global variables stream specified by 'output_filename'
 * and "PICCO_Global.txt". It then sets up various data structures for processing
 * the abstract syntax tree 'tree'. After processing the tree, it closes the global variables stream
 * and reads the content of "PICCO_Global.txt" into memory. Finally, it closes the output stream and
 * inserts the content of "PICCO_Global.txt" into 'output_filename' after a specific line.
 *
 * @param tree The abstract syntax tree to display.
 * @param output_filename The name of the output file where the generated code will be displayed.
 */
char *ast_show(aststmt tree, char *output_filename) {
    // This is where the output stream gets set 
    output = fopen(output_filename, "w+");
    // This is where the global variables stream get set 
    global_stream = fopen("PICCO_Global.txt", "w+");
    if (global_stream == NULL) {
        perror("Error opening/creating PICCO_Global.txt");
        return 1;
    }
    _curfile = Symbol(filename);
    indlev = 0;
    if_top = if_stack_new();
    if_tree = if_branchtree_new();
    bcs = batch_condition_stack_new();
    bss = batch_statement_stack_new();
    bpis = batch_private_index_stack_new();
    sns = struct_node_stack_new();
    ast_stmt_show(tree, if_tree);
    fclose(global_stream);
    fclose(output);

    // Read the global_stream as a string 
    char *fileContent = readFileToString();
    if_branchtree_free(if_tree);
    // Return the global_stream to picco.c
    return fileContent;
}

void ast_expr_show_stderr(astexpr tree) {
    int backup2 = dup(2); /* duplicate stderr */
    dup2(2, 1);           /* stdout to stderr */

    if (tree)
        ast_expr_show(tree);

    close(2);
    dup2(backup2, 2); /* restore stderr */
    close(backup2);
}

void ast_spec_show_stderr(astspec tree) {
    int backup = dup(1); /* duplicate stdout */
    dup2(2, 1);          /* stderr to stdout */

    if (tree)
        ast_spec_show(tree);

    dup2(backup, 1); /* restore stdout */
    close(backup);
}

void ast_decl_show_stderr(astdecl tree) {
    int backup = dup(1); /* duplicate stdout */
    dup2(2, 1);          /* stderr to stdout */

    if (tree)
        ast_decl_show(tree);

    dup2(backup, 1); /* restore stdout */
    close(backup);
}

void ast_show_stderr(aststmt tree) {
    int backup = dup(1); /* duplicate stdout */
    dup2(2, 1);          /* stderr to stdout */

    // ast_show(tree);

    dup2(backup, 1); /* restore stdout */
    close(backup);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                               *
 *     OMPi-EXTENSION NODES                                      *
 *                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void ast_oxclause_show(oxclause t) {
    if (t->type == OX_OCLIST) {
        if (t->u.list.next != NULL) {
            ast_oxclause_show(t->u.list.next);
            fprintf(output, " ");
        }
        assert((t = t->u.list.elem) != NULL);
    }

    fprintf(output, "%s", oxclausenames[t->type]);
    switch (t->type) {
    case OX_OCREDUCE:
           fprintf(output, "(%s : ", clausesubs[t->operator]);
           ast_decl_show(t->u.varlist);
           fprintf(output, ")");
           break;
    case OX_OCIN:
    case OX_OCOUT:
    case OX_OCINOUT:
           fprintf(output, "(");
           ast_decl_show(t->u.varlist);
           fprintf(output, ")");
           break;
    case OX_OCATNODE:
    case OX_OCATWORKER:
    case OX_OCSTART:
    case OX_OCSTRIDE:
           fprintf(output, "(");
           ast_expr_show(t->u.expr);
           fprintf(output, ")");
           break;
    case OX_OCSCOPE:
           fprintf(output, "scope(%s)", t->u.value == OX_SCOPE_NODES ? "nodes" : t->u.value == OX_SCOPE_WLOCAL ? "workers,local"
                                                                             : t->u.value == OX_SCOPE_WGLOBAL  ? "workers,global"
                                                                                                               : "???");
           break;
    }
}

void ast_oxdir_show(oxdir t) {
    fprintf(output, "#pragma ompix %s ", oxdirnames[t->type]);
    if (t->clauses)
           ast_oxclause_show(t->clauses);
    fprintf(output, "\n");
}

void ast_oxcon_show(oxcon t, branchnode current) {
    ast_oxdir_show(t->directive);
    if (t->body) {
           indent();
           ast_stmt_show(t->body, current);
    }
}
