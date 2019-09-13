
import numpy as np

START_FROM = [
    'random',
    'form_one',
    'uniform'
]

ANT_SYS_MODIF = [
    'cycle',
    'density',
    'quantity'
]


class AntAlgorithm:

    def __init__(self,
                 a=1,  # коэффициент запаха
                 b=1,  # коэффициент расстояния
                 rho=0.5,  # коэффициент высыхания
                 Q=1,  # количество выпускаемого феромона
                 e=0,  # количество элитных муравьев
                 local_refresh=False,  # локальное обновление феромона
                 where_to_start='random',  # место начала движения муравьев
                 as_modif='density',  # модификация обновления феромона
                 random_seed=-1  # фиксация генератора случайных чисел
                 ):
        self.a = a
        self.b = b
        self.rho = rho
        self.Q = Q
        self.e = int(e)
        self.local_refresh = local_refresh
        if where_to_start in START_FROM:
            self.where_to_start = where_to_start
        else:
            raise ValueError('wrong where to start')
        if as_modif in ANT_SYS_MODIF:
            self.as_modif = as_modif
        else:
            raise ValueError('wrong ant-system modification name')
        if random_seed > 0:
            np.random.seed(random_seed)

    # START_FROM = {'random': randomly,
            #   'form_one': from_one,
            #   'uniform': uniform}

    def fit(self,
            L,  # Матрица смежности графа
            AGES=-1,  # количество поколений
            ANTS=-1,  # количество муравьев в поколении
            ph=-1  # начальное значение ферромона

            ):
        self.L = L
        self.CITIES = len(L)
        if AGES > 0:
            self.AGES = int(AGES)
        else:
            self.AGES = self.CITIES * 50
        if ANTS > 0:
            self.ANTS = int(ANTS)
        else:
            self.ANTS = self.CITIES
        if ph >= 0:
            self.ph = ph
        else:
            self.ph = self.Q/self.CITIES

        # инициализация матрицы "краткости" дуг графа
        rev_L = 1/L
        # инициализация матрицы феромонов
        tao = np.ones((self.CITIES, self.CITIES)) * self.ph

        self.BEST_DIST = float("inf")                # лучшая длина маршрута
        self.BEST_ROUTE = None                       # лучший маршрут
        # матрица маршрутов муравьев в одном поколении (номера узлов графа)
        antROUTE = np.zeros((self.ANTS, self.CITIES))
        # вектор длины маршрута муравьев в одном поколении
        antDIST = np.zeros(self.ANTS)
        # вектор лучших длин маршрутов в каждом поколении
        self.antBEST_DIST = np.zeros(self.AGES)
        self.antAVERAGE_DIST = np.zeros(self.AGES)

        # основной цикл алгоритма
        # ---------- начало освновного цикла ----------
        for age in range(self.AGES):
            antROUTE.fill(0)
            antDIST.fill(0)

            # ---------- начало цикла обхода графа муравьями ----------
            for k in range(self.ANTS):

                if self.where_to_start == 'random':
                    # начальное расположение муравья в графе (случайное)
                    antROUTE[k, 0] = np.random.randint(
                        low=0, high=self.CITIES-1, size=1)
                elif self.where_to_start == 'uniform':
                    # начальное расположение муравья в графе (равномерное)
                    antROUTE[k, 0] = k % self.CITIES
                elif self.where_to_start == 'from_one':
                    # начальное расположение муравья в графе (все с одного)
                    antROUTE[k, 0] = 1
                else:
                    assert 'wrong where to start'

                # ---------- начало обхода графа k-ым муравьем ----------
                for s in range(1, self.CITIES):
                    # текущее положение муравья
                    from_city = int(antROUTE[k, s-1])
                    P = (tao[from_city] ** self.a) * \
                        (rev_L[from_city] ** self.b)
                    # вероятность посещения уже посещенных городов = 0
                    for i in range(s):
                        P[int(antROUTE[k, i])] = 0

                    # вероятность выбора направления, сумма всех P = 1
                    assert (np.sum(P) > 0), \
                        "Division by zero. P = %s,"\
                        " \n tao = %s \n rev_L = %s" % (
                        P, tao[from_city], rev_L[from_city])
                    P = P / np.sum(P)
                    # выбираем направление
                    isNotChosen = True
                    while isNotChosen:
                        rand = np.random.rand()
                        for p, to in zip(P, list(range(self.CITIES))):
                            if p >= rand:
                                # записываем город №s в вектор k-ого муравья
                                antROUTE[k, s] = to
                                isNotChosen = False
                                break
                    # локальное обновление феромона
                    if self.local_refresh:
                        for s in range(self.CITIES):
                            city_to = int(antROUTE[k, s])
                            city_from = int(antROUTE[k, s-1])
                            if self.as_modif == 'cycle':
                                tao[city_from, city_to] = \
                                    tao[city_from, city_to] + \
                                    (self.Q / antDIST[k])
                            elif self.as_modif == 'density':
                                tao[city_from, city_to] = tao[city_from,
                                                              city_to] + self.Q
                            elif self.as_modif == 'quantity':
                                tao[city_from, city_to] = \
                                    tao[city_from, city_to] + \
                                    (self.Q / L[city_from, city_to])
                            else:
                                assert 'wrong ant-system modification name'
                            tao[city_to, city_from] = tao[city_from, city_to]
                # ---------- конец цила обхода графа ----------

                # вычисляем длину маршрута k-ого муравья
                for i in range(self.CITIES):
                    city_from = int(antROUTE[k, i-1])
                    city_to = int(antROUTE[k, i])
                    antDIST[k] += L[city_from, city_to]

                # сравниваем длину маршрута с лучшим показателем
                if antDIST[k] < self.BEST_DIST:
                    self.BEST_DIST = antDIST[k]
                    self.BEST_ROUTE = antROUTE[k]
            # ---------- конец цикла обхода графа муравьями ----------

            # ---------- обновление феромонов----------
            # высыхание по всем маршрутам (дугам графа)
            tao *= (1-self.rho)

            # цикл обновления феромона
            for k in range(self.ANTS):
                for s in range(self.CITIES):
                    city_to = int(antROUTE[k, s])
                    city_from = int(antROUTE[k, s-1])
                    if self.as_modif == 'cycle':
                        tao[city_from, city_to] = \
                            tao[city_from, city_to] + \
                            (self.Q / antDIST[k])
                    elif self.as_modif == 'density':
                        tao[city_from, city_to] = \
                            tao[city_from, city_to] + self.Q
                    elif self.as_modif == 'quantity':
                        tao[city_from, city_to] = \
                            tao[city_from, city_to] + \
                            (self.Q / L[city_from, city_to])
                    else:
                        assert 'wrong ant-system modification name'
                    tao[city_to, city_from] = tao[city_from, city_to]

            # проход элитных е-муравьев по лучшему маршруту
            if self.e > 0:
                for s in range(self.CITIES):
                    city_to = int(self.BEST_ROUTE[s])
                    city_from = int(self.BEST_ROUTE[s-1])
                    if self.as_modif == 'cycle':
                        tao[city_from, city_to] = \
                            tao[city_from, city_to] + \
                            ((self.Q * self.e) / self.BEST_DIST)
                    elif self.as_modif == 'density':
                        tao[city_from, city_to] = \
                            tao[city_from, city_to] + self.Q * self.e
                    elif self.as_modif == 'quantity':
                        tao[city_from, city_to] = \
                            tao[city_from, city_to] + \
                            ((self.Q * self.e) / L[city_from, city_to])
                    else:
                        assert 'wrong ant-system modification name'
                    tao[city_to, city_from] = tao[city_from, city_to]

            # ---------- конец обновления феромона ----------

            # конец поколения муравьев

            # сбор информации для графиков
            self.antBEST_DIST[age] = self.BEST_DIST
            self.antAVERAGE_DIST[age] = np.average(antDIST)

        return tao

    def get_best_route(self):
        return self.BEST_ROUTE

    def get_best_dist(self):
        return self.BEST_DIST

    def get_best_dists(self):
        return [age for age in range(self.AGES)], self.antBEST_DIST

    def get_average_dists(self):
        return [age for age in range(self.AGES)], self.antAVERAGE_DIST