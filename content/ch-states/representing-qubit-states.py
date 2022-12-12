#!/usr/bin/env python
# coding: utf-8

# # Representing Qubit States

# You now know something about bits, and about how our familiar digital computers work. All the complex variables, objects and data structures used in modern software are basically all just big piles of bits. Those of us who work on quantum computing call these *classical variables.* The computers that use them, like the one you are using to read this article, we call *classical computers*.
# 
# In quantum computers, our basic variable is the _qubit:_ a quantum variant of the bit. These have exactly the same restrictions as normal bits do: they can store only a single binary piece of information, and can only ever give us an output of `0` or `1`. However, they can also be manipulated in ways that can only be described by quantum mechanics. This gives us new gates to play with, allowing us to find new ways to design algorithms.
# 
# To fully understand these new gates, we first need to understand how to write down qubit states. For this we will use the mathematics of vectors, matrices, and complex numbers. Though we will introduce these concepts as we go, it would be best if you are comfortable with them already. If you need a more in-depth explanation or a refresher, you can find the guide [here](https://qiskit.org/textbook/ch-appendix/linear_algebra.html).
# 
# 
# 
# 
# ## Contents
# 
# 1. [Classical vs Quantum Bits](#cvsq)    
#     1.1 [Statevectors](#statevectors)      
#     1.2 [Qubit Notation](#notation)    
#     1.3 [Exploring Qubits with Qiskit](#exploring-qubits)    
# 2. [The Rules of Measurement](#rules-measurement)    
#     2.1 [A Very Important Rule](#important-rule)    
#     2.2 [The Implications of this Rule](#implications)
# 3. [The Bloch Sphere](#bloch-sphere)    
#     3.1 [Describing the Restricted Qubit State](#bloch-sphere-1)     
#     3.2 [Visually Representing a Qubit State](#bloch-sphere-2)     
# 

# ## 1. Classical vs Quantum Bits <a id="cvsq"></a>
# 
# ### 1.1 Statevectors<a id="statevectors"></a>
# 
# In quantum physics we use _statevectors_ to describe the state of our system. Say we wanted to describe the position of a car along a track, this is a classical system so we could use a number $x$:
# 
# ![tracking a car with scalars](images/car_track_1.jpg)
# 
# $$ x=4 $$
# 
# Alternatively, we could instead use a collection of numbers in a vector called a  _statevector._ Each element in the statevector contains the probability of finding the car in a certain place:
# 
# ![tracking a car with vectors](images/car_track_2.jpg)
# 
# $$
# |x\rangle = \begin{bmatrix} 0\\ \vdots \\ 0 \\ 1 \\ 0 \\ \vdots \\ 0 \end{bmatrix} 
#             \begin{matrix} \\  \\  \\ \leftarrow \\  \\  \\  \\ \end{matrix}
#              \begin{matrix} \\  \\ \text{Probability of} \\ \text{car being at} \\ \text{position 4} \\  \\  \\ \end{matrix}   
# $$
# 
# This isn’t limited to position, we could also keep a statevector of all the possible speeds the car could have, and all the possible colours the car could be. With classical systems (like the car example above), this is a silly thing to do as it requires keeping huge vectors when we only really need one number. But as we will see in this chapter, statevectors happen to be a very good way of keeping track of quantum systems, including quantum computers.
# 
# 
# ### 1.2 Qubit Notation <a id="notation"></a>
# 
# Classical bits always have a completely well-defined state: they are either `0` or `1` at every point during a computation. There is no more detail we can add to the state of a bit than this. So to write down the state of a of classical bit (`c`), we can just use these two binary values. For example:
# 
#     c = 0
# 
# This restriction is lifted for quantum bits. Whether we get a `0` or a `1` from a qubit only needs to be well-defined when a measurement is made to extract an output. At that point, it must commit to one of these two options. At all other times, its state will be something more complex than can be captured by a simple binary value.
# 
# To see how to describe these, we can first focus on the two simplest cases. As we saw in the last section, it is possible to prepare a qubit in a state for which it definitely gives the outcome `0` when measured.
# 
# We need a name for this state. Let's be unimaginative and call it $0$ . Similarly, there exists a qubit state that is certain to output a `1`. We'll call this $1$. These two states are completely mutually exclusive. Either the qubit definitely outputs a ```0```, or it definitely outputs a ```1```. There is no overlap. One way to represent this with mathematics is to use two orthogonal vectors.
# 
# $$
# |0\rangle = \begin{bmatrix} 1 \\ 0 \end{bmatrix} \, \, \, \, |1\rangle =\begin{bmatrix} 0 \\ 1 \end{bmatrix}.
# $$
# 
# This is a lot of notation to take in all at once. First, let's unpack the weird $|$ and $\rangle$. Their job is essentially just to remind us that we are talking about the vectors that represent qubit states labelled $0$ and $1$. This helps us distinguish them from things like the bit values ```0``` and ```1``` or the numbers 0 and 1. It is part of the bra-ket notation, introduced by Dirac.
# 
# If you are not familiar with vectors, you can essentially just think of them as lists of numbers which we manipulate using certain rules. If you are familiar with vectors from your high school physics classes, you'll know that these rules make vectors well-suited for describing quantities with a magnitude and a direction. For example, the velocity of an object is described perfectly with a vector. However, the way we use vectors for quantum states is slightly different to this, so don't hold on too hard to your previous intuition. It's time to do something new!
# 
# With vectors we can describe more complex states than just $|0\rangle$ and $|1\rangle$. For example, consider the vector
# 
# $$
# |q_0\rangle = \begin{bmatrix} \tfrac{1}{\sqrt{2}} \\ \tfrac{i}{\sqrt{2}} \end{bmatrix} .
# $$
# 
# To understand what this state means, we'll need to use the mathematical rules for manipulating vectors. Specifically, we'll need to understand how to add vectors together and how to multiply them by scalars.
# 
# <p>
#  <details>
#   <summary>Reminder: Matrix Addition and Multiplication by Scalars (Click here to expand)</summary>
#   <p>To add two vectors, we add their elements together:
#     $$|a\rangle = \begin{bmatrix}a_0 \\ a_1 \\ \vdots \\ a_n \end{bmatrix}, \quad
#     |b\rangle = \begin{bmatrix}b_0 \\ b_1 \\ \vdots \\ b_n \end{bmatrix}$$
#     $$|a\rangle + |b\rangle = \begin{bmatrix}a_0 + b_0 \\ a_1 + b_1 \\ \vdots \\ a_n + b_n \end{bmatrix} $$
#     </p>
#   <p>And to multiply a vector by a scalar, we multiply each element by the scalar:
#     $$x|a\rangle = \begin{bmatrix}x \times a_0 \\ x \times  a_1 \\ \vdots \\ x \times  a_n \end{bmatrix}$$
#     </p>
#   <p>These two rules are used to rewrite the vector $|q_0\rangle$ (as shown above):
#       $$
#       \begin{aligned} 
#       |q_0\rangle & = \tfrac{1}{\sqrt{2}}|0\rangle + \tfrac{i}{\sqrt{2}}|1\rangle \\
#                   & = \tfrac{1}{\sqrt{2}}\begin{bmatrix}1\\0\end{bmatrix} + \tfrac{i}{\sqrt{2}}\begin{bmatrix}0\\1\end{bmatrix}\\
#                   & = \begin{bmatrix}\tfrac{1}{\sqrt{2}}\\0\end{bmatrix} + \begin{bmatrix}0\\\tfrac{i}{\sqrt{2}}\end{bmatrix}\\
#                   & = \begin{bmatrix}\tfrac{1}{\sqrt{2}} \\ \tfrac{i}{\sqrt{2}} \end{bmatrix}\\
#       \end{aligned}
#       $$
#  </details>
# </p>
# <p>
#  <details>
#   <summary>Reminder: Orthonormal Bases (Click here to expand)</summary>
#   <p>
#       It was stated before that the two vectors $|0\rangle$ and $|1\rangle$ are orthonormal, this means they are both <i>orthogonal</i> and <i>normalised</i>. Orthogonal means the vectors are at right angles:
#   </p><p><img src="images/basis.svg"></p>
#   <p>And normalised means their magnitudes (length of the arrow) is equal to 1. The two vectors $|0\rangle$ and $|1\rangle$ are <i>linearly independent</i>, which means we cannot describe $|0\rangle$ in terms of $|1\rangle$, and vice versa. However, using both the vectors $|0\rangle$ and $|1\rangle$, and our rules of addition and multiplication by scalars, we can describe all possible vectors in 2D space:
#     </p><p><img src="images/basis2.svg"></p>
#   <p>Because the vectors $|0\rangle$ and $|1\rangle$ are linearly independent, and can be used to describe any vector in 2D space using vector addition and scalar multiplication, we say the vectors $|0\rangle$ and $|1\rangle$ form a <i>basis</i>. In this case, since they are both orthogonal and normalised, we call it an <i>orthonormal basis</i>.
#  </details>
# </p>
# 
# Since the states $|0\rangle$ and $|1\rangle$ form an orthonormal basis, we can represent any 2D vector with a combination of these two states. This allows us to write the state of our qubit in the alternative form:
# 
# $$ |q_0\rangle = \tfrac{1}{\sqrt{2}}|0\rangle + \tfrac{i}{\sqrt{2}}|1\rangle $$
# 
# This vector, $|q_0\rangle$ is called the qubit's _statevector,_ it tells us everything we could possibly know about this qubit. For now, we are only able to draw a few simple conclusions about this particular example of a statevector: it is not entirely $|0\rangle$ and not entirely $|1\rangle$. Instead, it is described by a linear combination of the two. In quantum mechanics, we typically describe linear combinations such as this using the word 'superposition'.
# 
# Though our example state $|q_0\rangle$ can be expressed as a superposition of $|0\rangle$ and $|1\rangle$, it is no less a definite and well-defined qubit state than they are. To see this, we can begin to explore how a qubit can be manipulated.
# 
# ### 1.3 Exploring Qubits with Qiskit <a id="exploring-qubits"></a>
# 
# First, we need to import all the tools we will need:

# In[1]:


from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram, plot_bloch_vector
from math import sqrt, pi


# In Qiskit, we use the `QuantumCircuit` object to store our circuits, this is essentially a list of the quantum operations on our circuit and the qubits they are applied to.

# In[2]:


qc = QuantumCircuit(1) # Create a quantum circuit with one qubit


# In our quantum circuits, our qubits always start out in the state $|0\rangle$. We can use the `initialize()` method to transform this into any state. We give `initialize()` the vector we want in the form of a list, and tell it which qubit(s) we want to initialize in this state:

# In[3]:


qc = QuantumCircuit(1)  # Create a quantum circuit with one qubit
initial_state = [0,1]   # Define initial_state as |1>
qc.initialize(initial_state, 0) # Apply initialisation operation to the 0th qubit
qc.draw()  # Let's view our circuit


# We can then use one of Qiskit’s simulators to view the resulting state of our qubit.

# In[4]:


sim = Aer.get_backend('aer_simulator')  # Tell Qiskit how to simulate our circuit


# To get the results from our circuit, we use `run` to execute our circuit, giving the circuit and the backend as arguments. We then use `.result()` to get the result of this:

# In[5]:


qc = QuantumCircuit(1)  # Create a quantum circuit with one qubit
initial_state = [0,1]   # Define initial_state as |1>
qc.initialize(initial_state, 0) # Apply initialisation operation to the 0th qubit
qc.save_statevector()   # Tell simulator to save statevector
qobj = assemble(qc)     # Create a Qobj from the circuit for the simulator to run
result = sim.run(qobj).result() # Do the simulation and return the result


# from `result`, we can then get the final statevector using `.get_statevector()`:

# In[6]:


out_state = result.get_statevector()
print(out_state) # Display the output state vector


# **Note:** Python uses `j` to represent $i$ in complex numbers. We see a vector with two complex elements: `0.+0.j` = 0, and `1.+0.j` = 1.
# 
# Let’s now measure our qubit as we would in a real quantum computer and see the result:

# In[7]:


qc.measure_all()
qc.draw()


# This time, instead of the statevector we will get the counts for the `0` and `1` results using `.get_counts()`:

# In[8]:


qobj = assemble(qc)
result = sim.run(qobj).result()
counts = result.get_counts()
plot_histogram(counts)


# We can see that we (unsurprisingly) have a 100% chance of measuring $|1\rangle$. This time, let’s instead put our qubit into a superposition and see what happens. We will use the state $|q_0\rangle$ from earlier in this section:
# 
# $$ |q_0\rangle = \tfrac{1}{\sqrt{2}}|0\rangle + \tfrac{i}{\sqrt{2}}|1\rangle $$
# 
# We need to add these amplitudes to a python list. To add a complex amplitude, Python uses `j` for the imaginary unit (we normally call it "$i$" mathematically):

# In[9]:


initial_state = [1/sqrt(2), 1j/sqrt(2)]  # Define state |q_0>


# And we then repeat the steps for initialising the qubit as before:

# In[10]:


qc = QuantumCircuit(1) # Must redefine qc
qc.initialize(initial_state, 0) # Initialize the 0th qubit in the state `initial_state`
qc.save_statevector() # Save statevector
qobj = assemble(qc)
state = sim.run(qobj).result().get_statevector() # Execute the circuit
print(state)           # Print the result


# In[11]:


qobj = assemble(qc)
results = sim.run(qobj).result().get_counts()
plot_histogram(results)


# We can see we have equal probability of measuring either $|0\rangle$ or $|1\rangle$. To explain this, we need to talk about measurement.
# 
# ## 2. The Rules of Measurement <a id="rules-measurement"></a>
# ### 2.1 A Very Important Rule <a id="important-rule"></a>
# 
# There is a simple rule for measurement. To find the probability of measuring a state $|\psi \rangle$ in the state $|x\rangle$ we do:
# 
# $$p(|x\rangle) = | \langle x| \psi \rangle|^2$$
# 
# The symbols $\langle$ and $|$ tell us $\langle x |$ is a row vector and the symbols $|$ and $\rangle$ tell us $|\psi\rangle$ is a column vector. In quantum mechanics we call the column vectors _kets_ and the row vectors _bras._ Together they make up _bra-ket_ notation. Any ket $|a\rangle$ has a corresponding bra $\langle a|$, and we convert between them using the conjugate transpose.
# <details>
#     <summary>Reminder: Conjugate Transpose (Click here to expand)</summary>
#     <p>Conversion between bra-ket takes places using the <i>conjugate transpose</i> method. We know a ket (column vector) is represented as follows:
#         $$\quad|a\rangle = \begin{bmatrix}a_0 \\ a_1 \\ \vdots \\ a_n \end{bmatrix}$$
#     </p>
#     <p>
# In conjugate transpose method, the matrix is transposed and the elements are complex conjugated (represented by the "$*$" operation) where complex conjugate ("$*$") of a complex number is a number with an equal real part and an imaginary part equal in magnitude but opposite in sign. This gives the coressponding bra (row vector) as follows:
#         $$\langle a| = \begin{bmatrix}a_0^*, & a_1^*, & \dots & a_n^* \end{bmatrix}$$
#     </p>
# </details>   
# 
# <details>
#   <summary>Reminder: The Inner Product (Click here to expand)</summary>
#     <p>There are different ways to multiply vectors, here we use the <i>inner product</i>. The inner product is a generalisation of the <i>dot product</i> which you may already be familiar with. In this guide, we use the inner product between a bra (row vector) and a ket (column vector), and it follows this rule:
#         
# $$\langle a| = \begin{bmatrix}a_0^*, & a_1^*, & \dots & a_n^* \end{bmatrix}, \quad
#     |b\rangle = \begin{bmatrix}b_0 \\ b_1 \\ \vdots \\ b_n \end{bmatrix}$$
#     $$\langle a|b\rangle = a_0^* b_0 + a_1^* b_1 \dots a_n^* b_n$$
#     </p>
#   <p>We can see that the inner product of two vectors always gives us a scalar. A useful thing to remember is that the inner product of two orthogonal vectors is 0, for example if we have the orthogonal vectors $|0\rangle$ and $|1\rangle$:
#     $$\langle1|0\rangle = \begin{bmatrix} 0 & 1\end{bmatrix}\begin{bmatrix}1 \\ 0\end{bmatrix} = 0$$
#     </p>
#   <p>Additionally, remember that the vectors $|0\rangle$ and $|1\rangle$ are also normalised (magnitudes are equal to 1):
#     
# $$
#       \begin{aligned} 
#       \langle0|0\rangle & = \begin{bmatrix} 1 & 0\end{bmatrix}\begin{bmatrix}1 \\ 0\end{bmatrix} = 1 \\
#       \langle1|1\rangle & = \begin{bmatrix} 0 & 1\end{bmatrix}\begin{bmatrix}0 \\ 1\end{bmatrix} = 1
#       \end{aligned}
# $$
#    </p>
# </details>
# 
#          
# In the equation above, $|x\rangle$ can be any qubit state. To find the probability of measuring $|x\rangle$, we take the inner product of $|x\rangle$ and the state we are measuring (in this case $|\psi\rangle$), then square the magnitude. This may seem a little convoluted, but it will soon become second nature.
# 
# If we look at the state $|q_0\rangle$ from before, we can see the probability of measuring $|0\rangle$ is indeed $0.5$:
# 
# $$
# \begin{aligned}
# |q_0\rangle & = \tfrac{1}{\sqrt{2}}|0\rangle + \tfrac{i}{\sqrt{2}}|1\rangle \\
# \langle 0| q_0 \rangle  & = \tfrac{1}{\sqrt{2}}\langle 0|0\rangle + \tfrac{i}{\sqrt{2}}\langle 0|1\rangle \\
# & = \tfrac{1}{\sqrt{2}}\cdot 1 +  \tfrac{i}{\sqrt{2}} \cdot 0\\
# & = \tfrac{1}{\sqrt{2}}\\
# |\langle 0| q_0 \rangle|^2 & = \tfrac{1}{2}
# \end{aligned}
# $$
# 
# You should verify the probability of measuring $|1\rangle$ as an exercise.
# 
# This rule governs how we get information out of quantum states. It is therefore very important for everything we do in quantum computation. It also immediately implies several important facts.
# 
# ### 2.2 The Implications of this Rule <a id="implications"></a>
# ### #1 Normalisation
# 
# The rule shows us that amplitudes are related to probabilities. If we want the probabilities to add up to 1 (which they should!), we need to ensure that the statevector is properly normalized. Specifically, we need the magnitude of the state vector to be 1.
# 
# $$ \langle\psi|\psi\rangle = 1 \\ $$
# 
# Thus if:
# 
# $$ |\psi\rangle = \alpha|0\rangle + \beta|1\rangle $$
# 
# Then:
# 
# $$ |\alpha|^2 + |\beta|^2 = 1 $$
# 
# This explains the factors of $\sqrt{2}$ you have seen throughout this chapter. In fact, if we try to give `initialize()` a vector that isn’t normalised, it will give us an error:

# In[12]:


vector = [1,1]
# qc.initialize(vector, 0)


# #### Quick Exercise
# 1. Create a state vector that will give a $1/3$ probability of measuring $|0\rangle$.
# 2. Create a different state vector that will give the same measurement probabilities.
# 3. Verify that the probability of measuring $|1\rangle$ for these two states is $2/3$.

# You can check your answer in the widget below (accepts answers ±1% accuracy, you can use numpy terms such as '`pi`' and '`sqrt()`' in the vector):

# In[13]:


# Run the code in this cell to interact with the widget
from qiskit_textbook.widgets import state_vector_exercise
state_vector_exercise(target=1/3)


# ### #2 Alternative measurement
# 
# The measurement rule gives us the probability $p(|x\rangle)$ that a state $|\psi\rangle$ is measured as $|x\rangle$. Nowhere does it tell us that $|x\rangle$ can only be either $|0\rangle$ or $|1\rangle$.
# 
# The measurements we have considered so far are in fact only one of an infinite number of possible ways to measure a qubit. For any orthogonal pair of states, we can define a measurement that would cause a qubit to choose between the two.
# 
# This possibility will be explored more in the next section. For now, just bear in mind that $|x\rangle$ is not limited to being simply $|0\rangle$ or $|1\rangle$.

# ### #3 Global Phase
# 
# We know that measuring the state $|1\rangle$ will give us the output `1` with certainty. But we are also able to write down states such as 
# 
# $$\begin{bmatrix}0 \\ i\end{bmatrix} = i|1\rangle.$$
# 
# To see how this behaves, we apply the measurement rule.
# 
# $$ |\langle x| (i|1\rangle) |^2 = | i \langle x|1\rangle|^2 = |\langle x|1\rangle|^2 $$
# 
# Here we find that the factor of $i$ disappears once we take the magnitude of the complex number. This effect is completely independent of the measured state $|x\rangle$. It does not matter what measurement we are considering, the probabilities for the state $i|1\rangle$ are identical to those for $|1\rangle$. Since measurements are the only way we can extract any information from a qubit, this implies that these two states are equivalent in all ways that are physically relevant.
# 
# More generally, we refer to any overall factor $\gamma$ on a state for which $|\gamma|=1$ as a 'global phase'. States that differ only by a global phase are physically indistinguishable.
# 
# $$ |\langle x| ( \gamma |a\rangle) |^2 = | \gamma \langle x|a\rangle|^2 = |\langle x|a\rangle|^2 $$
# 
# Note that this is distinct from the phase difference _between_ terms in a superposition, which is known as the 'relative phase'. This becomes relevant once we consider different types of measurement and multiple qubits.
# 
# 
# ### #4 The Observer Effect
# 
# We know that the amplitudes contain information about the probability of us finding the qubit in a specific state, but once we have measured the qubit, we know with certainty what the state of the qubit is. For example, if we measure a qubit in the state:
# 
# $$ |q\rangle = \alpha|0\rangle + \beta|1\rangle$$
# 
# And find it in the state $|0\rangle$, if we measure again, there is a 100% chance of finding the qubit in the state $|0\rangle$. This means the act of measuring _changes_ the state of our qubits.
# 
# $$ |q\rangle = \begin{bmatrix} \alpha \\ \beta \end{bmatrix} \xrightarrow{\text{Measure }|0\rangle} |q\rangle = |0\rangle = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$$
# 
# We sometimes refer to this as _collapsing_ the state of the qubit. It is a potent effect, and so one that must be used wisely. For example, were we to constantly measure each of our qubits to keep track of their value at each point in a computation, they would always simply be in a well-defined state of either $|0\rangle$ or $|1\rangle$. As such, they would be no different from classical bits and our computation could be easily replaced by a classical computation. To achieve truly quantum computation we must allow the qubits to explore more complex states. Measurements are therefore only used when we need to extract an output. This means that we often place all the measurements at the end of our quantum circuit. 
# 
# We can demonstrate this using Qiskit’s statevector simulator. Let's initialize a qubit in superposition:

# In[14]:


qc = QuantumCircuit(1) # We are redefining qc
initial_state = [0.+1.j/sqrt(2),1/sqrt(2)+0.j]
qc.initialize(initial_state, 0)
qc.draw()


# This should initialize our qubit in the state:
# 
# $$ |q\rangle = \tfrac{i}{\sqrt{2}}|0\rangle + \tfrac{1}{\sqrt{2}}|1\rangle $$
# 
# We can verify this using the simulator:

# In[15]:


qc.save_statevector()
result = sim.run(assemble(qc)).result()
state = result.get_statevector()
print("Qubit State = " + str(state))


# We can see here the qubit is initialized in the state `[0.+0.70710678j 0.70710678+0.j]`, which is the state we expected.
# 
# Let’s now create a circuit where we measure this qubit:

# In[16]:


qc = QuantumCircuit(1) # We are redefining qc
initial_state = [0.+1.j/sqrt(2),1/sqrt(2)+0.j]
qc.initialize(initial_state, 0)
qc.measure_all()
qc.save_statevector()
qc.draw()


# When we simulate this entire circuit, we can see that one of the amplitudes is _always_ 0:

# In[17]:


qobj = assemble(qc)
state = sim.run(qobj).result().get_statevector()
print("State of Measured Qubit = " + str(state))


# You can re-run this cell a few times to reinitialize the qubit and measure it again. You will notice that either outcome is equally probable, but that the state of the qubit is never a superposition of $|0\rangle$ and $|1\rangle$. Somewhat interestingly, the global phase on the state $|0\rangle$ survives, but since this is global phase, we can never measure it on a real quantum computer.
# 
# ### A Note about Quantum Simulators
# 
# We can see that writing down a qubit’s state requires keeping track of two complex numbers, but when using a real quantum computer we will only ever receive a yes-or-no (`0` or `1`) answer for each qubit. The output of a 10-qubit quantum computer will look like this:
# 
# `0110111110`
# 
# Just 10 bits, no superposition or complex amplitudes. When using a real quantum computer, we cannot see the states of our qubits mid-computation, as this would destroy them! This behaviour is not ideal for learning, so Qiskit provides different quantum simulators: By default, the `aer_simulator` mimics the execution of a real quantum computer, but will also allow you to peek at quantum states before measurement if we include certain instructions in our circuit. For example, here we have included the instruction `.save_statevector()`, which means we can use `.get_statevector()` on the result of the simulation. 
# 
# 
# 

# ## 3. The Bloch Sphere <a id="bloch-sphere"></a>
# ### 3.1 Describing the Restricted Qubit State <a id="bloch-sphere-1"></a>
# 
# We saw earlier in this chapter that the general state of a qubit ($|q\rangle$) is:
# 
# $$
# |q\rangle = \alpha|0\rangle + \beta|1\rangle
# $$
# 
# $$
# \alpha, \beta \in \mathbb{C}
# $$
# 
# (The second line tells us $\alpha$ and $\beta$ are complex numbers). The first two implications in section 2 tell us that we cannot differentiate between some of these states. This means we can be more specific in our description of the qubit. 
# 
# Firstly, since we cannot measure global phase, we can only measure the difference in phase between the states $|0\rangle$ and $|1\rangle$. Instead of having $\alpha$ and $\beta$ be complex, we can confine them to the real numbers and add a term to tell us the relative phase between them:
# 
# $$
# |q\rangle = \alpha|0\rangle + e^{i\phi}\beta|1\rangle
# $$
# 
# $$
# \alpha, \beta, \phi \in \mathbb{R}
# $$
# 
# Finally, since the qubit state must be normalised, i.e.
# 
# $$
# \sqrt{\alpha^2 + \beta^2} = 1
# $$
# 
# we can use the trigonometric identity:
# 
# $$
# \sqrt{\sin^2{x} + \cos^2{x}} = 1
# $$
# 
# to describe the real $\alpha$ and $\beta$ in terms of one variable, $\theta$:
# 
# $$
# \alpha = \cos{\tfrac{\theta}{2}}, \quad \beta=\sin{\tfrac{\theta}{2}}
# $$
# 
# From this we can describe the state of any qubit using the two variables $\phi$ and $\theta$:
# 
# $$
# |q\rangle = \cos{\tfrac{\theta}{2}}|0\rangle + e^{i\phi}\sin{\tfrac{\theta}{2}}|1\rangle
# $$
# 
# $$
# \theta, \phi \in \mathbb{R}
# $$
# 
# ### 3.2 Visually Representing a Qubit State <a id="bloch-sphere-2"></a>
# 
# We want to plot our general qubit state:
# 
# $$
# |q\rangle = \cos{\tfrac{\theta}{2}}|0\rangle + e^{i\phi}\sin{\tfrac{\theta}{2}}|1\rangle
# $$
# 
# If we interpret $\theta$ and $\phi$ as spherical co-ordinates ($r = 1$, since the magnitude of the qubit state is $1$), we can plot any single qubit state on the surface of a sphere, known as the _Bloch sphere._
# 
# Below we have plotted a qubit in the state $|{+}\rangle$. In this case, $\theta = \pi/2$ and $\phi = 0$.
# 
# (Qiskit has a function to plot a bloch sphere, `plot_bloch_vector()`, but at the time of writing it only takes cartesian coordinates. We have included a function that does the conversion automatically).
# 
# 
# You can also try [this interactive Bloch sphere demo](https://javafxpert.github.io/grok-bloch/).
# 

# In[19]:


from qiskit_textbook.widgets import plot_bloch_vector_spherical
coords = [pi/2,0,1] # [Theta, Phi, Radius]
plot_bloch_vector_spherical(coords) # Bloch Vector with spherical coordinates


# #### Warning!
# When first learning about qubit states, it's easy to confuse the qubits _statevector_ with its _Bloch vector_. Remember the statevector is the vector discussed in [1.1](#notation), that holds the amplitudes for the two states our qubit can be in. The Bloch vector is a visualisation tool that maps the 2D, complex statevector onto real, 3D space.

# #### Quick Exercise
# Use `plot_bloch_vector()` or `plot_bloch_vector_spherical()` to plot a qubit in the states:
# 1. $|0\rangle$
# 2. $|1\rangle$
# 3. $\tfrac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$
# 4. $\tfrac{1}{\sqrt{2}}(|0\rangle - i|1\rangle)$
# 5. $\tfrac{1}{\sqrt{2}}\begin{bmatrix}i\\1\end{bmatrix}$

# In[24]:


def plot_statevector(α, β) -> None:
    from math import acos
    θ = 2 * acos(abs(α))
    from cmath import phase
    φ = phase(β) - phase(α)
    print(θ, φ)
    from IPython.display import display
    display(plot_bloch_vector_spherical([θ, φ, 1]))
plot_statevector(1, 0)
plot_statevector(0, 1)
plot_statevector(1 / sqrt(2), 1 / sqrt(2))
plot_statevector(1 / sqrt(2), - 1j / sqrt(2))
plot_statevector(1j / sqrt(2), 1 / sqrt(2))


# We have also included below a widget that converts from spherical co-ordinates to cartesian, for use with `plot_bloch_vector()`:

# In[18]:


from qiskit_textbook.widgets import bloch_calc
bloch_calc()


# In[19]:


import qiskit.tools.jupyter
get_ipython().run_line_magic('qiskit_version_table', '')

