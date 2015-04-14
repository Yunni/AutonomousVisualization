import itertools

from plot import ScatterPlot


unaries = ["log", "sqrt", "exp"]
binaries = ["+", "-", "*", "/", "**"]


def unary_expr(feature, mapping):
    '''Generate all expressions with unary operations

    :param feature: The feature name
    :param expr: A boolean variables indicate whether to map the feature
    :return:
    '''
    yield "v['" + feature + "']"
    for unary in unaries:
        if mapping:
            yield unary + "(" + feature + ")"
        else:
            yield unary + "(v['" + feature + "'])"


class ExprSearch:
    '''A class that greedily search all the scatter plots.

    '''

    def __init__(self, features, limit):
        '''

        :param features: Name of the features
        :param limit: The upper bound of variables in both expressions
        :return:
        '''
        self.features = features
        self.count = 0
        self.best = ScatterPlot("", "")
        self.limit = limit

    def test_expr(self, x_expr, y_expr, data):
        '''Test if the expression pair has better score than our current pair.
        If yes, set our current best pair to this pair.

        :param x_expr: A string expression of x
        :param y_expr: A string expression of y
        :param data: A list of list containing the whole numerical dataset
        '''
        plot = ScatterPlot(x_expr, y_expr)
        if not plot.set_points(self.features, data):
            return
        score = plot.calculate_score2()
        self.count += 1
        if score > self.best.score:
            self.best = plot
            print(self.best.score)
        print("Expr " + str(self.count) + ":" + str((x_expr, y_expr)) + "\t" + str(score))
        print("\tCurrent best" + str(self.best.expr) + "\t" + str(self.best.score))

    def pair_search(self, data):
        '''Search for a pair of expressions which both x, y contains 1 variable

        :param data:A list of list containing the whole numerical dataset
        :return:
        '''
        for x, y in itertools.combinations(range(len(self.features)), 2):
            for x_expr, y_expr in itertools.product(unary_expr(self.features[x], False),
                                                    unary_expr(self.features[y], False)):
                self.test_expr(x_expr, y_expr, data)

    def additional_search(self, data):
        '''Given the best expression. Try to add another

        :param data: A list of list containing the whole numerical dataset
        :return:
        '''
        if not self.best.expr[0]:
            self.pair_search(data)
        x_expr, y_expr = self.best.expr
        for feature in self.features:
            if feature in x_expr or feature in y_expr:
                continue
            for expr in unary_expr(feature, False):
                for binary in binaries:
                    self.test_expr(x_expr + binary + expr, y_expr, data)
                    self.test_expr(expr + binary + x_expr, y_expr, data)
                    self.test_expr(x_expr, y_expr + binary + expr, data)
                    self.test_expr(x_expr, expr + binary + y_expr, data)

    def search(self, data):
        '''Search for the best expressions for x and y.

        :param data: A list of list containing the whole numerical dataset
        :return:
        '''
        self.pair_search(data)
        for i in range(self.limit - 2):
            self.additional_search(data)
        return self.best
