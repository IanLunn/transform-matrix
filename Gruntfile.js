module.exports = function(grunt) {

    grunt.initConfig({
        pkg: grunt.file.readJSON('package.json'),

        sass: {
          dist: {
            options: {
              style: 'expanded'
            },
            files: {
              'css/styles.css': 'scss/styles.scss'
            }
          }
        },


        autoprefixer: {
          options: {
            browsers: ['last 2 version', 'ie 7', 'ie 8', 'ie 9']
          },
          multiple_files: {
            expand: true,
            flatten: true,
            src: 'css/*.css',
            dest: 'css/'
          }
        },


        cssmin: {
          combine: {
            files: {
              'css/styles-min.css': ['css/styles.css']
            }
          }
        },


        concat: {
          dist: {
            src: [
              'js/partials/third-party.js',
              'js/init.js',
            ],
            dest: 'js/scripts.js',
          },
        },


        uglify: {
          my_target: {
            files: {
              'js/scripts-min.js':
                ['js/scripts.js']
            }
          }
        },


        watch: {
          options: {
            livereload: true,
          },

          scripts: {
            files: ['js/**/*.js'],
            tasks: ['concat', 'uglify'],
            options: {
              spawn: false
            }
          },

          files: {
            files: ['*.html'],
            options: {
              spawn: false
            }
          },

          css: {
            files: ['scss/styles.scss', 'css/styles.css'],
            tasks: ['sass', 'autoprefixer', 'cssmin'],
            options: {
              spawn: false
            }
          }
        },


        connect: {
          server: {
            options: {
              port: 8000,
              base: './'
            }
          }
        }

    });

    require('load-grunt-tasks')(grunt);

    grunt.registerTask('default', ['connect', 'watch']); //watch for local development
    grunt.registerTask('run', ['sass', 'autoprefixer', 'cssmin', 'concat', 'uglify']); //manual compile
};