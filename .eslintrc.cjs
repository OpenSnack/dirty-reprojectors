module.exports = {
    root: true,
    parser: '@typescript-eslint/parser',
    parserOptions: {
        ecmaVersion: 2020
    },
    plugins: ['@typescript-eslint'],
    extends: [
        'eslint:recommended',
        'plugin:@typescript-eslint/recommended'
    ],
    settings: {
        'import/resolver': {
            typescript: {},
        },
    },
    env: {
        node: true
    },
    rules: {
        '@typescript-eslint/indent': ['error', 4]
    }
};
