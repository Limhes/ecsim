FROM quay.io/pypa/manylinux2014_x86_64

COPY . /code

WORKDIR /code/python

RUN yum install -y eigen3 \
    && for PYBIN in /opt/python/*/bin; do "${PYBIN}/python" setup.py bdist_wheel; done \
    && for whl in dist/*.whl; do auditwheel repair "${whl}"; done \
    && mv /code/python/wheelhouse /build-output

