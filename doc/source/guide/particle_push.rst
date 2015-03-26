
==================
 Particle Pushing
==================

We use Boris Method :cite:`boris70`. This method is well described in
:cite:`birdsall85` (ch. 4-3) (or more recent reprint) which is a standrart in
Particl-in-Cell technique for its good numerical properties:

Pros
----

* Second order accurate
* Error is **not** cummulative (as compared to Runge-Kutta)
* Fields are evelautaed only once
  
Cons
----

* The method requies to process particle array twice: This is required when
  used in hybrid codes as we perform only half time-step advance.


Derivation
----------

We are solving equation:

.. math::
   \frac{\mathbf{v}^{n+1/2} - \mathbf{v}^{n-1/2}}{\triangle t} =
   \frac{q}{m} \left[ \mathbf{E}^n + \mathbf{v}^{n} \times \mathbf{B}^n\right]

where :math:`\mathbf{v}^{n}` can be expressed as an average of two times

.. math::
   \mathbf{v}^{n} = \frac{1}{2}(\mathbf{v}^{n+1/2} + \mathbf{v}^{n-1/2}).

Lets perform following substitutions:

.. math::
   \mathbf{v}^{n+1/2} = \mathbf{v}^{+} + \frac{q \mathbf{E}^n}{m}
   \frac{\triangle t}{2}

   \mathbf{v}^{n-1/2} = \mathbf{v}^{-} - \frac{q \mathbf{E}^n}{m}
   \frac{\triangle t}{2}

by doing so, we split the particle push into rotation around magnetic field
:math:`\mathbf{B}^n` and acceleration by electric field
:math:`\mathbf{E}^n`. The rotation step is governed by the equation:

.. math::
   \mathbf{v}^{+} - \mathbf{v}^{-} =
   \frac{\triangle t}{2} \frac{q}{m}
   (\mathbf{v}^{+} + \mathbf{v}^{-}) \times \mathbf{B}


Steps
-----

#. :math:`\mathbf{v}^- = \mathbf{v}^{n-1/2} + \frac{q \mathbf{E}^n}{m}
   \frac{\triangle t}{2}`

#. :math:`\mathbf{t} = \mathbf{B}^n \frac{q}{m} \frac{\triangle t}{2}`

#. :math:`\mathbf{s} = 2 \mathbf{t}/ (1 + t^2)` (from the requirement
   :math:`|\mathbf{v}^-|^2=|\mathbf{v}^+|^2`)

#. :math:`\mathbf{v}^\prime = \mathbf{v}^- + \mathbf{v}^- \times \mathbf{t}`

#. :math:`\mathbf{v}^+ = \mathbf{v}^- + \mathbf{v}^\prime \times \mathbf{s}`

#. :math:`\mathbf{v}^{n + 1/2} = \mathbf{v}^+ + \frac{q \mathbf{E}^n}{m}
   \frac{\triangle t}{2}`


References
----------

.. bibliography:: refs.bib
