import copy
import numpy as np
import pandas as pd
import random
from itertools import combinations
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from math import pi

# region ReadOnly
nbPositions = 52
nbCourse = 100
# endregion

# region Paramètres "secondaires"
# Taux d'abeilles créer à chaque génération
birthRate = 20  # math.floor(.25 * nbCourse)

# Taux de parent
nbParent = 20  # math.floor(.25 * nbCourse)

# Nb d'abeilles qui vont muter et qui ne seront pas remplacer
nbMutants = nbCourse - (nbParent / 2)
# endregion

# region Paramètres primaires
nbGeneration = 75000

# Nb de permutation de gène pour créer des mutants
mutationRate = 2

# Taille de la séquence de gène à permuter
mutationEffect = 1
# endregion


class AlgoGen:

    def __init__(self):
        global birthRate, nbParent, nbMutants
        # Position de la ruche (départ et arrivée)
        self.startPos = (500, 500)

        is_new_swarm = True if input("Démarrer à partir d'un nouvel essain ? (O/N)").upper() == "O" else False
        if is_new_swarm:
            # Données initiales qui correspondent à la position des points (des fleurs) dans un ordre X
            self.data = Utils.get_data(is_new_swarm, "bees")

            # Obtention de 99 autres abeilles random et de leurs parcours
            self.courses = self.get_random_courses()
        else:
            self.courses = Utils.get_data(is_new_swarm, "last_output")
            self.courses = self.courses.reshape((nbCourse, -1, 2))

        # On ajuste le nombre de parent pour simplifer la suite du programme qui fonctionne par paire
        nbParent = 2 if nbParent == 1 or nbParent == 0 else nbParent - 1 if nbParent % 2 != 0 else nbParent
        nbMutants = nbCourse - birthRate - nbParent
        # Génération
        self.generate()

        Utils.plot_overview()

        if input("Sauvegarde ?(O/N) ").upper() == "O":
            # Sauvegarde
            Utils.set_data(self.courses)
            print("Sauvegardé")

    # Génère aléatoirement n parcours selon un parcours initial
    # Retour sous d'une liste np sous forme [[parcours1], [parcours2], ...]]
    def get_random_courses(self):
        courses = []
        for i in range(nbCourse):
            new_list = copy.deepcopy(self.data)
            # On conserve les données initial qui peuvent correspondre à des données déjà travaillées
            if i != 0:
                np.random.shuffle(new_list)

            new_list = np.concatenate([[self.startPos], new_list, [self.startPos]])
            courses.append(new_list)
        return np.array(courses)

    # Génération
    def generate(self):
        # region Initialization
        global mutationEffect, mutationRate
        n = -1
        step = nbGeneration / 100
        nb_slack, max_slack = 0, 50
        mutation_rates = random.choices([x for x in range(2, 5)], [0.6, 0.25, 0.15], k=10 ** 5)
        mutation_effects = random.choices([x for x in range(4, 8)], [0.25, 0.25, 0.25, 0.25], k=10 ** 5)
        distances = np.array([0 for _ in range(nbCourse)])
        # distances totale pour chaque course
        # []
        childs = [[[0 for _ in range(2)] for _ in range(nbPositions)] for _ in range(int(birthRate))]

        distances = self.fitness(distances)
        # endregion
        while n < nbGeneration:
            # region Debug
            n += 1
            percentage = n / step
            if percentage % 10 == 0:
                print("\nProgress:", str(percentage) + "%")
                print(mutationRate, mutationEffect)
                Utils.update_stats(self.courses, distances, n)
            last_min = distances[0]
            # Reset lorsqu'on suprasse la stagnation: retour aux vleurs par défaut (mR, mE)
            if last_min - distances[0] > 250:
                mutationRate = 2
                mutationEffect = 1
                nb_slack = 0
            if n % max_slack == 0:
                if last_min - distances[0] < 250:
                    nb_slack += 1
                    if nb_slack > max_slack:
                        mutationRate = 2
                        mutationEffect = 1
                        nb_slack = 0
                    else:
                        if nb_slack < max_slack * .25:
                            mutationEffect = 2
                        elif nb_slack < max_slack * .50:
                            mutationEffect = 3
                        else:
                            mutationEffect = mutation_effects.pop()
                        mutationRate = mutation_rates.pop()
            # endregion
            childs = self.crossover(childs)
            self.mutation()
            self.selection(childs)
            distances = self.fitness(distances)

    # Fonction d'évaluation qui calcul la distance entre chaque segment et la distance total pour une course
    # Retourne une liste np de même taille que la liste d'entrée sous la forme [d_seg1, d_seg2, ...]
    def fitness(self, distances):
        current_pos = self.startPos
        coef_man_to_pytha = 4/pi
        for j in range(nbCourse):
            total_distance = 0
            for i in range(nbPositions - 1):
                x, y = current_pos
                dest_x, dest_y = self.courses[j, i + 1, 0], self.courses[j, i + 1, 1]
                distance_segment = (abs(x - dest_x) + abs(y - dest_y))/coef_man_to_pytha
                total_distance += distance_segment
                current_pos = (dest_x, dest_y)
            distances[j] = total_distance

        self.courses = self.courses[distances.argsort()]
        return np.sort(distances)

    # Fonction qui retourne n enfants en provenance des p meilleurs abeilles
    # Les enfants sont génétiquement modifiés à partir des génomes parents
    def crossover(self, childs):
        parents_combinations = list(combinations([x for x in range(nbParent)], 2))
        np.random.shuffle(parents_combinations)
        split = int(nbPositions / 2)
        for i in range(birthRate):
            idx1, idx2 = parents_combinations[i]
            parents = [self.courses[idx1], self.courses[idx2]]
            genome_child = np.reshape([parents[0][:split], parents[1][split:]], (nbPositions, -1))
            for y, gene in enumerate(genome_child[1:split]):
                if gene in genome_child[split:-1]:
                    if random.random() > .5:
                        genitor = parents[0]
                    else:
                        genitor = parents[1]
                    for _, genitor_gen in enumerate(genitor[1:-1]):
                        if genitor_gen not in genome_child:
                            genome_child[y + 1] = genitor_gen
                            break
            childs[i] = genome_child

        return childs

    # Fonction qui gère la mutation appliqué aux abeilles:
    # Un nombre n de gènes est permuté pour obtenir de nouveaux individus
    # On tire 2 nombres au hasard (idx1, idx2) et on permute les séquences idx1->mutationEffect, idx2->mutationEffect
    def mutation(self):
        for i in range(int(nbParent / 2), nbCourse):
            mutated_course = copy.deepcopy(self.courses[i])
            for j in range(mutationRate):
                # region Déclarations d'index pour permuter les séquences
                idx1, idx2 = random.randrange(1, nbPositions - 1 - mutationEffect), \
                             random.randrange(1, nbPositions - 1 - mutationEffect)
                if idx1 == idx2:
                    if idx1 > 1:
                        idx1 -= 1
                    else:
                        idx2 += 1
                if idx1 > idx2:
                    idx1, idx2 = idx2, idx1
                # if idx2 + mutationEffect > nbPositions - 1:
                #     idx2 -= mutationEffect - (nbPositions - 1 - idx2)
                diff = idx2 - idx1
                if diff < mutationEffect:
                    if idx1 - (mutationEffect - diff) > 0:
                        idx1 -= mutationEffect - diff
                    elif idx2 + (mutationEffect - diff) < nbPositions - mutationEffect:
                        idx2 += mutationEffect - diff
                # endregion
                for w in range(mutationEffect):
                    mutated_course[[idx1 + w, idx2 + w]] = mutated_course[[idx2 + w, idx1 + w]]
            self.courses[i] = mutated_course

    # Fonction qui remplace les abeilles
    def selection(self, childs):
        for i in range(birthRate):
            self.courses[-(i + 1)] = childs[i]


class Utils:
    overview_min = np.zeros((0, 2), int)
    overview_avg = np.zeros((0, 2), int)

    @staticmethod
    def set_data(courses):
        df = pd.DataFrame({"x": np.ravel(courses[:, :, 0]), "y": np.ravel(courses[:, :, 1])})
        df.to_excel("data/output/last_output.xlsx", index=False)

    @staticmethod
    def get_data(is_new_swarm, fname):
        if is_new_swarm:
            return pd.read_excel("data/input/" + fname + ".xlsx").to_numpy()
        else:
            return pd.read_excel("data/output/" + fname + ".xlsx").to_numpy()

    @staticmethod
    def update_stats(courses, distances, percentage):
        min_distance = int(distances[0])
        avg_distance = int(np.mean(distances))

        Utils.update_overview(min_distance, avg_distance, percentage)
        Utils.plot_min_course(courses, min_distance, avg_distance)

    @staticmethod
    def plot_min_course(courses, min_distance, avg_distance):
        min_course = [courses[0][:, 0], courses[0][:, 1]]
        plt.plot(min_course[0], min_course[1], linewidth=2, linestyle=":", c="orange")
        plt.scatter(min_course[0], min_course[1], color="orange", linewidths=1, marker="o")

        plt.legend(["meilleur trajet: " + str(min_distance) + "m", "distance moyenne: " + str(avg_distance) + "m"],
                   markerscale=0,
                   frameon=False, borderaxespad=0, bbox_to_anchor=(0, 1.02, 1, .2), ncol=1, loc="lower left",
                   mode="expand")
        plt.scatter(500, 500, color="red", marker="p", s=200)
        plt.show()

    @staticmethod
    def plot_overview():
        x, y = Utils.overview_avg[:, 0], Utils.overview_avg[:, 1]
        plt.plot(x, y)
        x, y = Utils.overview_min[:, 0], Utils.overview_min[:, 1]
        plt.plot(x, y)
        orange_patch = mpatches.Patch(color='orange', label='Trajet minimum')
        blue_patch = mpatches.Patch(color='blue', label='Trajet moyen')

        plt.legend([orange_patch, blue_patch], ["trajet minimum", "trajet moyen"],
                   markerscale=0,
                   frameon=False, borderaxespad=0, bbox_to_anchor=(0, 1.02, 1, .2), ncol=2, loc="lower right")
        plt.suptitle('Evolution de la population au fil des générations', fontsize=16)
        plt.xlabel('N Générations', fontsize=12)
        plt.ylabel('Evolution', fontsize=12)

        plt.show()

    @staticmethod
    def update_overview(min_distance, avg_distance, percentage):
        Utils.overview_min = np.vstack((Utils.overview_min, [percentage, 1 / min_distance]))
        Utils.overview_avg = np.vstack((Utils.overview_avg, [percentage, 1 / avg_distance]))


AlgoGen()
